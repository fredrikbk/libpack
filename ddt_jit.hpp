#ifndef DDT_JIT_HPP
#define DDT_JIT_HPP

#include <map>
#include <string>
#include <vector>
#include <cstdio>

#include <mpi.h>

#define TRUNK 1
#ifdef TRUNK
#include "llvm/IR/Module.h"
#include "llvm/IR/LLVMContext.h"
#include "llvm/IR/DerivedTypes.h"
#include "llvm/IR/IRBuilder.h"
#include "llvm/IR/Intrinsics.h"
#else
#include "llvm/Module.h"
#include "llvm/LLVMContext.h"
#include "llvm/DerivedTypes.h"
#include "llvm/IRBuilder.h"
#include "llvm/Intrinsics.h"
#endif

#include "llvm/Support/TargetSelect.h"
#include "llvm/Analysis/Verifier.h"
#include "llvm/Analysis/Passes.h"
#include "llvm/ExecutionEngine/ExecutionEngine.h"
#include "llvm/ExecutionEngine/JIT.h"
#include "llvm/PassManager.h"
#include "llvm/LinkAllPasses.h"
#include "llvm/Transforms/Scalar.h"
#include "llvm/Transforms/IPO/PassManagerBuilder.h"

using namespace llvm;

static Function* func; // printf function for debugging
static Module *TheModule;
static IRBuilder<> Builder(getGlobalContext());
static std::map<std::string, Value*> NamedValues;
static ExecutionEngine *TheExecutionEngine;
//static FunctionPassManager *TheFPM;

/* Utility Functions */
Value* multNode(int op1, Value* op2PtrNode) {
    Value* op1Node = ConstantInt::get(getGlobalContext(), APInt(64, op1, false));
    Value* op2Node = Builder.CreateIntCast(op2PtrNode, Type::getInt64Ty(getGlobalContext()), false); 
    return Builder.CreateMul(op1Node, op2Node);
}

Value* constNode(int val) {
    return ConstantInt::get(getGlobalContext(), APInt(32, val, false));
}

/* Base class for all datatypes */
class FARC_Datatype {

    public:
    virtual ~FARC_Datatype() {}
    virtual Value *Codegen_Pack(Value* inbuf, Value* incount, Value* outbuf) = 0;
    virtual Value *Codegen_Unpack(Value* inbuf, Value* incount, Value* outbuf) = 0;
    virtual int getExtend() = 0;
    virtual int getSize() = 0;

    void* packer;
    void* unpacker;

};

/* Class for primitive types, such as MPI_INT, MPI_BYTE, etc */
class FARC_PrimitiveDatatype : public FARC_Datatype {

    MPI_Datatype Type; // this MUST be a primitive type
    int Extend;
    int Size;

    public:
    FARC_PrimitiveDatatype(MPI_Datatype type);
    Value* Codegen_Pack(Value* inbuf, Value* incount, Value* outbuf);
    Value* Codegen_Unpack(Value* inbuf, Value* incount, Value* outbuf);
    int getExtend();
    int getSize();

};

int FARC_PrimitiveDatatype::getExtend() {
    return this->Extend;
}

int FARC_PrimitiveDatatype::getSize() {
    return this->Size;
}

Value* FARC_PrimitiveDatatype::Codegen_Pack(Value* inbuf, Value* incount, Value* outbuf) {
                            
    /** for debugging */ 
    //Value *fmt_ptr1 = Builder.CreateGlobalStringPtr("in PTR: %p out PTR: %p size: %i\n\n\0");
    // now we can print as follows:
    //CallInst *call = Builder.CreateCall3(func, fmt_ptr1, inbuf, outbuf);

    Value* contig_extend = multNode(this->Extend, incount);
    Value* memcopy = Builder.CreateMemCpy(outbuf, inbuf, contig_extend, 1);

    // return 0 for now
    return Constant::getNullValue(Type::getInt32Ty(getGlobalContext()));

}

Value* FARC_PrimitiveDatatype::Codegen_Unpack(Value* inbuf_ptr, Value* incount, Value* outbuf_ptr) {

    // does exactly the same as pack for primitive types
    return this->Codegen_Pack(inbuf_ptr, incount, outbuf_ptr);
}

FARC_PrimitiveDatatype::FARC_PrimitiveDatatype(MPI_Datatype type) {

    this->Type = type;

    if (Type == MPI_BYTE)   this->Extend = 1;
    if (Type == MPI_CHAR)   this->Extend = 1;
    if (Type == MPI_DOUBLE) this->Extend = sizeof(double);
    if (Type == MPI_INT)    this->Extend = sizeof(int);
    //TODO add more

    this->Size = this->Extend;

}

/* Class for contiguous types */
class FARC_ContiguousDatatype : public FARC_Datatype {

    FARC_Datatype* Basetype;
    int Count;

    public:
    FARC_ContiguousDatatype(FARC_Datatype* type, int count) : Basetype(type), Count(count) {}
    int getExtend();
    int getSize();
    Value *Codegen_Pack(Value* inbuf, Value* incount, Value* outbuf);
    Value *Codegen_Unpack(Value* inbuf, Value* incount, Value* outbuf);

};

int FARC_ContiguousDatatype::getExtend() {
    return this->Count * this->Basetype->getExtend();
}

int FARC_ContiguousDatatype::getSize() {
    return this->Count * this->Basetype->getSize();
}

Value* FARC_ContiguousDatatype::Codegen_Pack(Value* inbuf, Value* incount, Value* outbuf) {

    Function* TheFunction = Builder.GetInsertBlock()->getParent();

    // Loop
    BasicBlock* PreheaderBB = Builder.GetInsertBlock();
    BasicBlock* LoopBB = BasicBlock::Create(getGlobalContext(), "loop", TheFunction);
    Builder.CreateBr(LoopBB);
    Builder.SetInsertPoint(LoopBB);

    // Induction var phi nodes
    PHINode *out = Builder.CreatePHI(Type::getInt8PtrTy(getGlobalContext()), 2, "out");
    out->addIncoming(outbuf, PreheaderBB);
    PHINode *in= Builder.CreatePHI(Type::getInt8PtrTy(getGlobalContext()), 2, "in");
    in->addIncoming(inbuf, PreheaderBB);
    PHINode *i = Builder.CreatePHI(Type::getInt32Ty(getGlobalContext()), 2, "i");
    i->addIncoming(constNode(0), PreheaderBB);


    // Basetype Code Generation
    Basetype->Codegen_Pack(in, ConstantInt::get(getGlobalContext(), APInt(32, this->Count, false)), out);


    // Increment the out ptr by Size(Basetype) * Blocklen
    Value* bytesWritten = ConstantInt::get(getGlobalContext(), APInt(64, this->Basetype->getSize() * this->Count, false));
    Value* out_addr_cvi = Builder.CreatePtrToInt(out, Type::getInt64Ty(getGlobalContext()));
    Value* out_addr = Builder.CreateAdd(out_addr_cvi, bytesWritten);
    Value* nextout = Builder.CreateIntToPtr(out_addr, Type::getInt8PtrTy(getGlobalContext()));

    // Increment the in ptr by Extend(Basetype) * Stride
    Value* bytesToStride = ConstantInt::get(getGlobalContext(), APInt(64, this->Basetype->getExtend() * this->Count, false));
    Value* in_addr_cvi = Builder.CreatePtrToInt(in, Type::getInt64Ty(getGlobalContext()));
    Value* in_addr = Builder.CreateAdd(in_addr_cvi, bytesToStride);
    Value* nextin = Builder.CreateIntToPtr(in_addr, Type::getInt8PtrTy(getGlobalContext()));

    // Increment outer loop index
    Value* nexti = Builder.CreateAdd(i, constNode(1), "nexti");
    Value* EndCond_outer = Builder.CreateICmpEQ(nexti, incount, "loopcond");

    // Create and branch to the outer loop postamble
    BasicBlock *LoopEndBB = Builder.GetInsertBlock();
    BasicBlock *AfterBB = BasicBlock::Create(getGlobalContext(), "afterloop", TheFunction);
    Builder.CreateCondBr(EndCond_outer, AfterBB, LoopBB);
    Builder.SetInsertPoint(AfterBB);

    // Add backedges for the outer loop induction variable
    out->addIncoming(nextout, LoopEndBB);
    in->addIncoming(nextin, LoopEndBB);
    i->addIncoming(nexti, LoopEndBB);

    return Constant::getNullValue(Type::getInt32Ty(getGlobalContext()));
}

Value* FARC_ContiguousDatatype::Codegen_Unpack(Value* inbuf, Value* incount, Value* outbuf) {

    Function* TheFunction = Builder.GetInsertBlock()->getParent();

    // Loop
    BasicBlock* PreheaderBB = Builder.GetInsertBlock();
    BasicBlock* LoopBB = BasicBlock::Create(getGlobalContext(), "loop", TheFunction);
    Builder.CreateBr(LoopBB);
    Builder.SetInsertPoint(LoopBB);

    // Induction var phi nodes
    PHINode *out = Builder.CreatePHI(Type::getInt8PtrTy(getGlobalContext()), 2, "out");
    out->addIncoming(outbuf, PreheaderBB);
    PHINode *in= Builder.CreatePHI(Type::getInt8PtrTy(getGlobalContext()), 2, "in");
    in->addIncoming(inbuf, PreheaderBB);
    PHINode *i = Builder.CreatePHI(Type::getInt32Ty(getGlobalContext()), 2, "i");
    i->addIncoming(constNode(0), PreheaderBB);


    // Basetype Code Generation
    Basetype->Codegen_Pack(in, ConstantInt::get(getGlobalContext(), APInt(32, this->Count, false)), out);


    // Increment the out ptr by Size(Basetype) * Blocklen
    Value* bytesWritten = ConstantInt::get(getGlobalContext(), APInt(64, this->Basetype->getExtend() * this->Count, false));
    Value* out_addr_cvi = Builder.CreatePtrToInt(out, Type::getInt64Ty(getGlobalContext()));
    Value* out_addr = Builder.CreateAdd(out_addr_cvi, bytesWritten);
    Value* nextout = Builder.CreateIntToPtr(out_addr, Type::getInt8PtrTy(getGlobalContext()));

    // Increment the in ptr by Extend(Basetype) * Stride
    Value* bytesToStride = ConstantInt::get(getGlobalContext(), APInt(64, this->Basetype->getSize() * this->Count, false));
    Value* in_addr_cvi = Builder.CreatePtrToInt(in, Type::getInt64Ty(getGlobalContext()));
    Value* in_addr = Builder.CreateAdd(in_addr_cvi, bytesToStride);
    Value* nextin = Builder.CreateIntToPtr(in_addr, Type::getInt8PtrTy(getGlobalContext()));

    // Increment outer loop index
    Value* nexti = Builder.CreateAdd(i, constNode(1), "nexti");
    Value* EndCond_outer = Builder.CreateICmpEQ(nexti, incount, "loopcond");

    // Create and branch to the outer loop postamble
    BasicBlock *LoopEndBB = Builder.GetInsertBlock();
    BasicBlock *AfterBB = BasicBlock::Create(getGlobalContext(), "afterloop", TheFunction);
    Builder.CreateCondBr(EndCond_outer, AfterBB, LoopBB);
    Builder.SetInsertPoint(AfterBB);

    // Add backedges for the outer loop induction variable
    out->addIncoming(nextout, LoopEndBB);
    in->addIncoming(nextin, LoopEndBB);
    i->addIncoming(nexti, LoopEndBB);

    return Constant::getNullValue(Type::getInt32Ty(getGlobalContext()));

}

/* Class for hvector types */
class FARC_HVectorDatatype : public FARC_Datatype {

    FARC_Datatype* Basetype;
    int Count;
    int Blocklen;
    int Stride;

    public:
    FARC_HVectorDatatype(FARC_Datatype* type, int count, int blocklen, int stride) : Basetype(type), Count(count), Blocklen(blocklen), Stride(stride) {}
    Value *Codegen_Pack(Value* inbuf, Value* incount, Value* outbuf);
    Value *Codegen_Unpack(Value* inbuf, Value* incount, Value* outbuf);
    int getExtend();
    int getSize();

};

int FARC_HVectorDatatype::getExtend() {

    return (this->Count-1)*this->Stride + this->Blocklen*this->Basetype->getExtend();

}

int FARC_HVectorDatatype::getSize() {

    return this->Count * this->Blocklen*this->Basetype->getSize();

}

Value* FARC_HVectorDatatype::Codegen_Pack(Value* inbuf_ptr, Value* incount, Value* outbuf_ptr) {

    /** for debugging **   
    std::vector<llvm::Type*> printf_arg_types;
    printf_arg_types.push_back(Type::getInt8PtrTy(getGlobalContext()));
    FunctionType* printf_type = FunctionType::get(Type::getInt32Ty(getGlobalContext()), printf_arg_types, true);
    Function *func = Function::Create(printf_type, Function::ExternalLinkage, Twine("printf"), TheModule);
    Value *fmt_ptr = Builder.CreateGlobalStringPtr("stride to add: %i\n\0");
    Value *fmt_ptr2 = Builder.CreateGlobalStringPtr("restore stride\n\0");
    // now we can print as follows:
    //llvm::CallInst *call = builder.CreateCall2(func, fmt_ptr, ValueToPrint);
    */
                            
    Value* inbuf = Builder.CreateLoad(inbuf_ptr);
    Value* outbuf = Builder.CreateLoad(outbuf_ptr);
                            
    // start value for the loop over incount
    Value *CounterValue_outer = ConstantInt::get(getGlobalContext(), APInt(32, 0, false));
                            
    // Make the new basic block for the loop header, inserting after current block.
    Function *TheFunction = Builder.GetInsertBlock()->getParent();
    BasicBlock *Preheader_outer_BB = Builder.GetInsertBlock();
    BasicBlock *Loop_outer_BB = BasicBlock::Create(getGlobalContext(), "loop", TheFunction);
                            
    // Insert an explicit fall through from the current block to the LoopBB.
    Builder.CreateBr(Loop_outer_BB);
                            
    // Start insertion in LoopBB.
    Builder.SetInsertPoint(Loop_outer_BB);
                            
    // Start the PHI node with an entry for Start.
    PHINode *Variable_outer = Builder.CreatePHI(Type::getInt32Ty(getGlobalContext()), 2, "i");
    Variable_outer->addIncoming(CounterValue_outer, Preheader_outer_BB);
                            
    // Emit the body of the outer loop.
                            
    Value *CounterValue_inner = ConstantInt::get(getGlobalContext(), APInt(32, 0, false));
                            
    // Make the new basic block for the loop header, inserting after current block.
    BasicBlock *Preheader_inner_BB = Builder.GetInsertBlock();
    BasicBlock *Loop_inner_BB = BasicBlock::Create(getGlobalContext(), "loop", TheFunction);
                            
    // Insert an explicit fall through from the current block to the LoopBB.
    Builder.CreateBr(Loop_inner_BB);
                            
    // Start insertion in LoopBB.
    Builder.SetInsertPoint(Loop_inner_BB);
                            
    // Start the PHI node with an entry for Start.
    PHINode *Variable_inner = Builder.CreatePHI(Type::getInt32Ty(getGlobalContext()), 2, "i");
    Variable_inner->addIncoming(CounterValue_inner, Preheader_inner_BB);
                            
    // emit the body for the inner loop
    // store oldin, so that we can calc stride
    Value* oldin = Builder.CreateLoad(inbuf_ptr);
                            
    Basetype->Codegen_Pack(inbuf_ptr, ConstantInt::get(getGlobalContext(), APInt(32, this->Blocklen, false)), outbuf_ptr);
                            
    Value* newin = Builder.CreateLoad(inbuf_ptr);
                            
    // *imbuf += stride - (newin - oldin)
    Value* oldin_int = Builder.CreatePtrToInt(oldin, Type::getInt64Ty(getGlobalContext()));
    Value* newin_int = Builder.CreatePtrToInt(newin, Type::getInt64Ty(getGlobalContext()));
    Value* written_bytes = Builder.CreateSub(newin_int, oldin_int);
    Value* stride = ConstantInt::get(getGlobalContext(), APInt(64, Stride, false));
    Value* st_add = Builder.CreateSub(stride, written_bytes);
                            
    Value* newin_with_stride_int = Builder.CreateAdd(st_add, newin_int);
    Value* newin_with_stride_ptr = Builder.CreateIntToPtr(newin_with_stride_int, Type::getInt8PtrTy(getGlobalContext()));
                            
    Builder.CreateStore(newin_with_stride_ptr, inbuf_ptr);
                            
    // Emit the step value. 
    Value* StepVal_inner = ConstantInt::get(getGlobalContext(), APInt(32, 1, false));
    Value* NextVar_inner = Builder.CreateAdd(Variable_inner, StepVal_inner, "nextvar");
                            
    // check if we are finished with the loop over count
    Value* count = ConstantInt::get(getGlobalContext(), APInt(32, Count, false));
    Value* EndCond_inner = Builder.CreateICmpNE(NextVar_inner, count, "loopcond");
                            
    // Create the "after loop" block and insert it.
    BasicBlock *LoopEnd_inner_BB = Builder.GetInsertBlock();
    BasicBlock *After_inner_BB = BasicBlock::Create(getGlobalContext(), "afterloop", TheFunction);
                            
    // Insert the conditional branch into the end of LoopEndBB.
    Builder.CreateCondBr(EndCond_inner, Loop_inner_BB, After_inner_BB);
                            
    // Any new code will be inserted in AfterBB.
    Builder.SetInsertPoint(After_inner_BB);
                            
    // Add a new entry to the PHI node for the backedge.
    Variable_inner->addIncoming(NextVar_inner, LoopEnd_inner_BB);
                            
    // Remove the stride that was added in the last iteration
    Builder.CreateStore(newin, inbuf_ptr);
                            
    /* llop epilog of outer loop*/
                            
    // Emit the step value. 
    Value* StepVal_outer = ConstantInt::get(getGlobalContext(), APInt(32, 1, false));
    Value* NextVar_outer = Builder.CreateAdd(Variable_outer, StepVal_outer, "nextvar");
    // check if we are finished with the loop over incount
    Value* EndCond_outer = Builder.CreateICmpNE(NextVar_outer, incount, "loopcond");
                            
    // Create the "after loop" block and insert it.
    BasicBlock *LoopEnd_outer_BB = Builder.GetInsertBlock();
    BasicBlock *After_outer_BB = BasicBlock::Create(getGlobalContext(), "afterloop", TheFunction);
    // Insert the conditional branch into the end of LoopEndBB.
    Builder.CreateCondBr(EndCond_outer, Loop_outer_BB, After_outer_BB);

    // Any new code will be inserted in AfterBB.
    Builder.SetInsertPoint(After_outer_BB);

    // Add a new entry to the PHI node for the backedge.
    Variable_outer->addIncoming(NextVar_outer, LoopEnd_outer_BB);

    // return 0 for now
    return Constant::getNullValue(Type::getInt32Ty(getGlobalContext()));


}

Value* FARC_HVectorDatatype::Codegen_Unpack(Value* inbuf_ptr, Value* incount, Value* outbuf_ptr) {

    /** for debugging *   
    std::vector<llvm::Type*> printf_arg_types;
    printf_arg_types.push_back(Type::getInt8PtrTy(getGlobalContext()));
    FunctionType* printf_type = FunctionType::get(Type::getInt32Ty(getGlobalContext()), printf_arg_types, true);
    Function *func = Function::Create(printf_type, Function::ExternalLinkage, Twine("printf"), TheModule);
    Value *fmt_ptr = Builder.CreateGlobalStringPtr("copy blocklen\n\0");
    // now we can print as follows:
    //llvm::CallInst *call = builder.CreateCall2(func, fmt_ptr, ValueToPrint);
    */
                            
    Value* inbuf = Builder.CreateLoad(inbuf_ptr);
    Value* outbuf = Builder.CreateLoad(outbuf_ptr);
                            
    // start value for the loop over incount
    Value *CounterValue_outer = ConstantInt::get(getGlobalContext(), APInt(32, 0, false));
                            
    // Make the new basic block for the loop header, inserting after current block.
    Function *TheFunction = Builder.GetInsertBlock()->getParent();
    BasicBlock *Preheader_outer_BB = Builder.GetInsertBlock();
    BasicBlock *Loop_outer_BB = BasicBlock::Create(getGlobalContext(), "loop", TheFunction);
                            
    // Insert an explicit fall through from the current block to the LoopBB.
    Builder.CreateBr(Loop_outer_BB);
                            
    // Start insertion in LoopBB.
    Builder.SetInsertPoint(Loop_outer_BB);
                            
    // Start the PHI node with an entry for Start.
    PHINode *Variable_outer = Builder.CreatePHI(Type::getInt32Ty(getGlobalContext()), 2, "i");
    Variable_outer->addIncoming(CounterValue_outer, Preheader_outer_BB);
                            
    // Emit the body of the outer loop.
                            
    Value *CounterValue_inner = ConstantInt::get(getGlobalContext(), APInt(32, 0, false));
                            
    // Make the new basic block for the loop header, inserting after current block.
    BasicBlock *Preheader_inner_BB = Builder.GetInsertBlock();
    BasicBlock *Loop_inner_BB = BasicBlock::Create(getGlobalContext(), "loop", TheFunction);
                            
    // Insert an explicit fall through from the current block to the LoopBB.
    Builder.CreateBr(Loop_inner_BB);
                            
    // Start insertion in LoopBB.
    Builder.SetInsertPoint(Loop_inner_BB);
                            
    // Start the PHI node with an entry for Start.
    PHINode *Variable_inner = Builder.CreatePHI(Type::getInt32Ty(getGlobalContext()), 2, "i");
    Variable_inner->addIncoming(CounterValue_inner, Preheader_inner_BB);
                            
    // emit the body for the inner loop
    // store oldout, so that we can calc stride
    Value* oldout = Builder.CreateLoad(outbuf_ptr);
                            
   // llvm::CallInst *call = Builder.CreateCall(func, fmt_ptr);
    Basetype->Codegen_Unpack(inbuf_ptr, ConstantInt::get(getGlobalContext(), APInt(32, this->Blocklen, false)), outbuf_ptr);
                            
    Value* newout = Builder.CreateLoad(outbuf_ptr);
                            
    // *outbuf += stride - (newout - oldout)
    Value* oldout_int = Builder.CreatePtrToInt(oldout, Type::getInt64Ty(getGlobalContext()));
    Value* newout_int = Builder.CreatePtrToInt(newout, Type::getInt64Ty(getGlobalContext()));
    Value* written_bytes = Builder.CreateSub(newout_int, oldout_int);
    Value* stride = ConstantInt::get(getGlobalContext(), APInt(64, this->Stride, false));
    Value* st_add = Builder.CreateSub(stride, written_bytes);
                            
    Value* newout_with_stride_int = Builder.CreateAdd(st_add, newout_int);
    Value* newout_with_stride_ptr = Builder.CreateIntToPtr(newout_with_stride_int, Type::getInt8PtrTy(getGlobalContext()));
                            
    Builder.CreateStore(newout_with_stride_ptr, outbuf_ptr);
                            
    // Emit the step value. 
    Value* StepVal_inner = ConstantInt::get(getGlobalContext(), APInt(32, 1, false));
    Value* NextVar_inner = Builder.CreateAdd(Variable_inner, StepVal_inner, "nextvar");
                            
    // check if we are finished with the loop over count
    Value* count = ConstantInt::get(getGlobalContext(), APInt(32, Count, false));
    Value* EndCond_inner = Builder.CreateICmpNE(NextVar_inner, count, "loopcond");
                            
    // Create the "after loop" block and insert it.
    BasicBlock *LoopEnd_inner_BB = Builder.GetInsertBlock();
    BasicBlock *After_inner_BB = BasicBlock::Create(getGlobalContext(), "afterloop", TheFunction);
                            
    // Insert the conditional branch into the end of LoopEndBB.
    Builder.CreateCondBr(EndCond_inner, Loop_inner_BB, After_inner_BB);
                            
    // Any new code will be inserted in AfterBB.
    Builder.SetInsertPoint(After_inner_BB);
                            
    // Add a new entry to the PHI node for the backedge.
    Variable_inner->addIncoming(NextVar_inner, LoopEnd_inner_BB);
                            
    // Remove the stride that was added in the last iteration
    Builder.CreateStore(newout, outbuf_ptr);
                            
    /* llop epilog of outer loop*/
                            
    // Emit the step value. 
    Value* StepVal_outer = ConstantInt::get(getGlobalContext(), APInt(32, 1, false));
    Value* NextVar_outer = Builder.CreateAdd(Variable_outer, StepVal_outer, "nextvar");
    // check if we are finished with the loop over incount
    Value* EndCond_outer = Builder.CreateICmpNE(NextVar_outer, incount, "loopcond");
                            
    // Create the "after loop" block and insert it.
    BasicBlock *LoopEnd_outer_BB = Builder.GetInsertBlock();
    BasicBlock *After_outer_BB = BasicBlock::Create(getGlobalContext(), "afterloop", TheFunction);
    // Insert the conditional branch into the end of LoopEndBB.
    Builder.CreateCondBr(EndCond_outer, Loop_outer_BB, After_outer_BB);

    // Any new code will be inserted in AfterBB.
    Builder.SetInsertPoint(After_outer_BB);

    // Add a new entry to the PHI node for the backedge.
    Variable_outer->addIncoming(NextVar_outer, LoopEnd_outer_BB);

    // return 0 for now
    return Constant::getNullValue(Type::getInt32Ty(getGlobalContext()));


}

/* Class for vector types */
class FARC_VectorDatatype : public FARC_Datatype {

    FARC_Datatype* Basetype;
    int Count;
    int Blocklen;
    int Stride;

    public:
    FARC_VectorDatatype(FARC_Datatype* type, int count, int blocklen, int stride) : Basetype(type), Count(count), Blocklen(blocklen), Stride(stride) {}
    Value *Codegen_Pack(Value* inbuf, Value* incount, Value* outbuf);
    Value *Codegen_Unpack(Value* inbuf, Value* incount, Value* outbuf);
    int getExtend();
    int getSize();

};

int FARC_VectorDatatype::getExtend() {

    return (this->Count - 1) * this->Basetype->getExtend() * this->Stride + this->Blocklen * this->Basetype->getExtend();

}

int FARC_VectorDatatype::getSize() {

    return this->Count * this->Blocklen*this->Basetype->getSize();

}

Value* FARC_VectorDatatype::Codegen_Pack(Value* inbuf, Value* incount, Value* outbuf) {

    /** for debugging **   
    std::vector<llvm::Type*> printf_arg_types;
    printf_arg_types.push_back(Type::getInt8PtrTy(getGlobalContext()));
    FunctionType* printf_type = FunctionType::get(Type::getInt32Ty(getGlobalContext()), printf_arg_types, true);
    Function *func = Function::Create(printf_type, Function::ExternalLinkage, Twine("printf"), TheModule);
    Value *fmt_ptr = Builder.CreateGlobalStringPtr("stride to add: %i\n\0");
    Value *fmt_ptr2 = Builder.CreateGlobalStringPtr("restore stride\n\0");
    // now we can print as follows:
    //llvm::CallInst *call = builder.CreateCall2(func, fmt_ptr, ValueToPrint);
    */

    Function *TheFunction = Builder.GetInsertBlock()->getParent();

    // Outer loop
    BasicBlock *Preheader_outer_BB = Builder.GetInsertBlock();
    BasicBlock *Loop_outer_BB = BasicBlock::Create(getGlobalContext(), "outerloop", TheFunction);
    Builder.CreateBr(Loop_outer_BB);
    Builder.SetInsertPoint(Loop_outer_BB);

    // Induction var phi nodes
    PHINode *out_outer = Builder.CreatePHI(Type::getInt8PtrTy(getGlobalContext()), 2, "out");
    out_outer->addIncoming(outbuf, Preheader_outer_BB);
    PHINode *in_outer= Builder.CreatePHI(Type::getInt8PtrTy(getGlobalContext()), 2, "in");
    in_outer->addIncoming(inbuf, Preheader_outer_BB);
    PHINode *i = Builder.CreatePHI(Type::getInt32Ty(getGlobalContext()), 2, "i");
    i->addIncoming(constNode(0), Preheader_outer_BB);
    

    
    // Inner loop
    BasicBlock *Preheader_inner_BB = Builder.GetInsertBlock();
    BasicBlock *Loop_inner_BB = BasicBlock::Create(getGlobalContext(), "innerloop", TheFunction);
    Builder.CreateBr(Loop_inner_BB);
    Builder.SetInsertPoint(Loop_inner_BB);

    // Induction var phi nodes
    PHINode *out_inner = Builder.CreatePHI(Type::getInt8PtrTy(getGlobalContext()), 2, "out1");
    out_inner->addIncoming(out_outer, Preheader_inner_BB);
    PHINode *in_inner= Builder.CreatePHI(Type::getInt8PtrTy(getGlobalContext()), 2, "in1");
    in_inner->addIncoming(in_outer, Preheader_inner_BB);
    PHINode *j = Builder.CreatePHI(Type::getInt32Ty(getGlobalContext()), 2, "j");
    j->addIncoming(constNode(0), Preheader_inner_BB);


    // Basetype Code Generation
    Basetype->Codegen_Pack(in_inner, ConstantInt::get(getGlobalContext(), APInt(32, this->Blocklen, false)), out_inner);


    // Increment the out ptr by Size(Basetype) * Blocklen
    Value* bytesWritten = ConstantInt::get(getGlobalContext(), APInt(64, this->Basetype->getSize() * this->Blocklen, false));
    Value* out_addr_cvi = Builder.CreatePtrToInt(out_inner, Type::getInt64Ty(getGlobalContext()));
    Value* out_addr = Builder.CreateAdd(out_addr_cvi, bytesWritten);
    Value* nextout_inner = Builder.CreateIntToPtr(out_addr, Type::getInt8PtrTy(getGlobalContext()));

    // Increment the in ptr by Extend(Basetype) * Stride
    Value* bytesToStride = ConstantInt::get(getGlobalContext(), APInt(64, this->Basetype->getExtend() * this->Stride, false));
    Value* in_addr_cvi = Builder.CreatePtrToInt(in_inner, Type::getInt64Ty(getGlobalContext()));
    Value* in_addr = Builder.CreateAdd(in_addr_cvi, bytesToStride);
    Value* nextin_inner = Builder.CreateIntToPtr(in_addr, Type::getInt8PtrTy(getGlobalContext()));

    // Prefetch
/*
    std::vector<Type *> arg_type;
    Function *fun = Intrinsic::getDeclaration(TheFunction->getParent(), Intrinsic::prefetch, arg_type);
    std::vector<Value*> args;
    args.push_back(inbuf);
    args.push_back(constNode(0));
    args.push_back(constNode(1));
    args.push_back(constNode(1));
    Builder.CreateCall(fun, args);
*/

    // Increment inner loop index
    Value* stepj = ConstantInt::get(getGlobalContext(), APInt(32, 1, false));
    Value* nextj = Builder.CreateAdd(j, constNode(1), "nextj");
    
    // check if we are finished with the loop over count
    Value* count = ConstantInt::get(getGlobalContext(), APInt(32, Count, false));
    Value* EndCond_inner = Builder.CreateICmpEQ(nextj, count, "innercond");
                            
    // Create and branch to the inner loop postamble
    BasicBlock *LoopEnd_inner_BB = Builder.GetInsertBlock();
    BasicBlock *After_inner_BB = BasicBlock::Create(getGlobalContext(), "afterinner", TheFunction);
    Builder.CreateCondBr(EndCond_inner, After_inner_BB, Loop_inner_BB);
    Builder.SetInsertPoint(After_inner_BB);

    // Add backedges for the inner loop induction variables
    out_inner->addIncoming(nextout_inner, LoopEnd_inner_BB);
    in_inner->addIncoming(nextin_inner, LoopEnd_inner_BB);
    j->addIncoming(nextj, LoopEnd_inner_BB);


    // Move the in ptr back Extend(Basetype) * Stride - Size(Basetype) * Blocklen  
    Value* nextin_outer_val = Builder.CreateAdd(in_addr_cvi, bytesWritten);
    Value* nextin_outer = Builder.CreateIntToPtr(nextin_outer_val, Type::getInt8PtrTy(getGlobalContext()));

    // Increment outer loop index
    Value* nexti = Builder.CreateAdd(i, constNode(1), "nexti");
    Value* EndCond_outer = Builder.CreateICmpEQ(nexti, incount, "loopcond");

    // Create and branch to the outer loop postamble
    BasicBlock *LoopEnd_outer_BB = Builder.GetInsertBlock();
    BasicBlock *After_outer_BB = BasicBlock::Create(getGlobalContext(), "afterloop", TheFunction);
    Builder.CreateCondBr(EndCond_outer, After_outer_BB, Loop_outer_BB);
    Builder.SetInsertPoint(After_outer_BB);

    // Add backedges for the outer loop induction variable
    out_outer->addIncoming(nextout_inner, LoopEnd_outer_BB);
    in_outer->addIncoming(nextin_outer, LoopEnd_outer_BB);
    i->addIncoming(nexti, LoopEnd_outer_BB);

    // return 0 for now
    return Constant::getNullValue(Type::getInt32Ty(getGlobalContext()));
}

Value* FARC_VectorDatatype::Codegen_Unpack(Value* inbuf_ptr, Value* incount, Value* outbuf_ptr) {

    /** for debugging **   
    std::vector<llvm::Type*> printf_arg_types;
    printf_arg_types.push_back(Type::getInt8PtrTy(getGlobalContext()));
    FunctionType* printf_type = FunctionType::get(Type::getInt32Ty(getGlobalContext()), printf_arg_types, true);
    Function *func = Function::Create(printf_type, Function::ExternalLinkage, Twine("printf"), TheModule);
    Value *fmt_ptr = Builder.CreateGlobalStringPtr("stride to add: %i\n\0");
    Value *fmt_ptr2 = Builder.CreateGlobalStringPtr("restore stride\n\0");
    // now we can print as follows:
    //llvm::CallInst *call = builder.CreateCall2(func, fmt_ptr, ValueToPrint);
    */
                            
    Value* inbuf = Builder.CreateLoad(inbuf_ptr);
    Value* outbuf = Builder.CreateLoad(outbuf_ptr);
                            
    // start value for the loop over incount
    Value *CounterValue_outer = ConstantInt::get(getGlobalContext(), APInt(32, 0, false));
                            
    // Make the new basic block for the loop header, inserting after current block.
    Function *TheFunction = Builder.GetInsertBlock()->getParent();
    BasicBlock *Preheader_outer_BB = Builder.GetInsertBlock();
    BasicBlock *Loop_outer_BB = BasicBlock::Create(getGlobalContext(), "loop", TheFunction);
                            
    // Insert an explicit fall through from the current block to the LoopBB.
    Builder.CreateBr(Loop_outer_BB);
                            
    // Start insertion in LoopBB.
    Builder.SetInsertPoint(Loop_outer_BB);
                            
    // Start the PHI node with an entry for Start.
    PHINode *Variable_outer = Builder.CreatePHI(Type::getInt32Ty(getGlobalContext()), 2, "i");
    Variable_outer->addIncoming(CounterValue_outer, Preheader_outer_BB);
                            
    // Emit the body of the outer loop.
                            
    Value *CounterValue_inner = ConstantInt::get(getGlobalContext(), APInt(32, 0, false));
                            
    // Make the new basic block for the loop header, inserting after current block.
    BasicBlock *Preheader_inner_BB = Builder.GetInsertBlock();
    BasicBlock *Loop_inner_BB = BasicBlock::Create(getGlobalContext(), "loop", TheFunction);
                            
    // Insert an explicit fall through from the current block to the LoopBB.
    Builder.CreateBr(Loop_inner_BB);
                            
    // Start insertion in LoopBB.
    Builder.SetInsertPoint(Loop_inner_BB);
                            
    // Start the PHI node with an entry for Start.
    PHINode *Variable_inner = Builder.CreatePHI(Type::getInt32Ty(getGlobalContext()), 2, "i");
    Variable_inner->addIncoming(CounterValue_inner, Preheader_inner_BB);
                            
    // emit the body for the inner loop
    // store oldout, so that we can calc stride
    Value* oldout = Builder.CreateLoad(outbuf_ptr);
                            
    Basetype->Codegen_Unpack(inbuf_ptr, ConstantInt::get(getGlobalContext(), APInt(32, this->Blocklen, false)), outbuf_ptr);
                            
    Value* newout = Builder.CreateLoad(outbuf_ptr);
                            
    // *outbuf += stride - (newout - oldout)
    Value* oldout_int = Builder.CreatePtrToInt(oldout, Type::getInt64Ty(getGlobalContext()));
    Value* newout_int = Builder.CreatePtrToInt(newout, Type::getInt64Ty(getGlobalContext()));
    Value* read_bytes = Builder.CreateSub(newout_int, oldout_int);
    Value* stride = ConstantInt::get(getGlobalContext(), APInt(64, Stride*this->Basetype->getExtend(), false));
    Value* st_add = Builder.CreateSub(stride, read_bytes);
                            
    Value* newout_with_stride_int = Builder.CreateAdd(st_add, newout_int);
    Value* newout_with_stride_ptr = Builder.CreateIntToPtr(newout_with_stride_int, Type::getInt8PtrTy(getGlobalContext()));
                            
    Builder.CreateStore(newout_with_stride_ptr, outbuf_ptr);
                            
    // Emit the step value. 
    Value* StepVal_inner = ConstantInt::get(getGlobalContext(), APInt(32, 1, false));
    Value* NextVar_inner = Builder.CreateAdd(Variable_inner, StepVal_inner, "nextvar");
                            
    // check if we are finished with the loop over count
    Value* count = ConstantInt::get(getGlobalContext(), APInt(32, Count, false));
    Value* EndCond_inner = Builder.CreateICmpNE(NextVar_inner, count, "loopcond");
                            
    // Create the "after loop" block and insert it.
    BasicBlock *LoopEnd_inner_BB = Builder.GetInsertBlock();
    BasicBlock *After_inner_BB = BasicBlock::Create(getGlobalContext(), "afterloop", TheFunction);
                            
    // Insert the conditional branch into the end of LoopEndBB.
    Builder.CreateCondBr(EndCond_inner, Loop_inner_BB, After_inner_BB);
                            
    // Any new code will be inserted in AfterBB.
    Builder.SetInsertPoint(After_inner_BB);
                            
    // Add a new entry to the PHI node for the backedge.
    Variable_inner->addIncoming(NextVar_inner, LoopEnd_inner_BB);
                            
    // Remove the stride that was added in the last iteration
    Builder.CreateStore(newout, outbuf_ptr);
                            
    /* llop epilog of outer loop*/
                            
    // Emit the step value. 
    Value* StepVal_outer = ConstantInt::get(getGlobalContext(), APInt(32, 1, false));
    Value* NextVar_outer = Builder.CreateAdd(Variable_outer, StepVal_outer, "nextvar");
    // check if we are finished with the loop over incount
    Value* EndCond_outer = Builder.CreateICmpNE(NextVar_outer, incount, "loopcond");
                            
    // Create the "after loop" block and insert it.
    BasicBlock *LoopEnd_outer_BB = Builder.GetInsertBlock();
    BasicBlock *After_outer_BB = BasicBlock::Create(getGlobalContext(), "afterloop", TheFunction);
    // Insert the conditional branch into the end of LoopEndBB.
    Builder.CreateCondBr(EndCond_outer, Loop_outer_BB, After_outer_BB);

    // Any new code will be inserted in AfterBB.
    Builder.SetInsertPoint(After_outer_BB);

    // Add a new entry to the PHI node for the backedge.
    Variable_outer->addIncoming(NextVar_outer, LoopEnd_outer_BB);

    // return 0 for now
    return Constant::getNullValue(Type::getInt32Ty(getGlobalContext()));


}

/* Class for struct types */
class FARC_StructDatatype : public FARC_Datatype {

    int Count;
    std::vector<FARC_Datatype*> Types;
    std::vector<int> Blocklen;
    std::vector<int> Displ;

    public:
    FARC_StructDatatype(int count, int* blocklen, long*  displ, FARC_Datatype** types);
    Value *Codegen_Pack(Value* inbuf, Value* incount, Value* outbuf);
    Value *Codegen_Unpack(Value* inbuf, Value* incount, Value* outbuf);
    int getExtend();
    int getSize();

};

FARC_StructDatatype::FARC_StructDatatype(int count, int* blocklen, long*  displ, FARC_Datatype** types) {

    this->Count = count;
    for (int i=0; i<count; i++) Blocklen.push_back(blocklen[i]);
    for (int i=0; i<count; i++) Displ.push_back(displ[i]);
    for (int i=0; i<count; i++) Types.push_back(types[i]);

}

int FARC_StructDatatype::getExtend() {

    if (this->Count == 0) return 0;

    int lb = this->Displ[0];
    int ub = this->Displ[0] + this->Types[0]->getExtend() * this->Blocklen[0];
    for (int i=0; i<this->Count; i++) {
        int tmp_ub = this->Displ[i] + this->Types[i]->getExtend() * this->Blocklen[i];
        int tmp_lb = this->Displ[i];
        if (tmp_ub > ub) ub = tmp_ub;
        if (tmp_lb < lb) lb = tmp_lb;
    }

    return ub-lb;

}

int FARC_StructDatatype::getSize() {

    int sum = 0;
    for (int i=0; i<this->Count; i++) {
        sum += this->Types[i]->getSize() * this->Blocklen[i];
    }

    return sum;

}

Value* FARC_StructDatatype::Codegen_Pack(Value* inbuf_ptr, Value* incount, Value* outbuf_ptr) {

    // start value for the loop over incount
    Value* CounterValue = ConstantInt::get(getGlobalContext(), APInt(32, 0, false));
    // save inbuf ptr so that we can calculate inbuf = inbuf_orig + i * extend for each iter over incount
    Value* inbuf_orig = Builder.CreateLoad(inbuf_ptr);
    Value* inbuf_orig_int = Builder.CreatePtrToInt(inbuf_orig, Type::getInt64Ty(getGlobalContext()));
    Value* extend = ConstantInt::get(getGlobalContext(), APInt(32, this->getExtend(), false));
                            
    // Make the new basic block for the loop header, inserting after current block.
    Function* TheFunction = Builder.GetInsertBlock()->getParent();
    BasicBlock* PreheaderBB = Builder.GetInsertBlock();
    BasicBlock* LoopBB = BasicBlock::Create(getGlobalContext(), "loop", TheFunction);
                            
    // Insert an explicit fall through from the current block to the LoopBB.
    Builder.CreateBr(LoopBB);
                            
    // Start insertion in LoopBB.
    Builder.SetInsertPoint(LoopBB);
                            
    // Start the PHI node with an entry for Start.
    PHINode* Variable = Builder.CreatePHI(Type::getInt32Ty(getGlobalContext()), 2, "i");
    Variable->addIncoming(CounterValue, PreheaderBB);
                            
    // Emit the body of the loop:

    //inbuf = inbuf_old + incount_i * extend
    Value* incnt_i_mul_ext_32 = Builder.CreateMul(Variable, extend);
    Value* incnt_i_mul_ext_64 = Builder.CreateZExt(incnt_i_mul_ext_32, Type::getInt64Ty(getGlobalContext()));
    Value* inbuf_new_int = Builder.CreateAdd(inbuf_orig_int, incnt_i_mul_ext_64);

    for (int i=0; i<this->Count; i++) {
        // inbuf_displ = inptr_new + displ[i]
        Value* displ_i = ConstantInt::get(getGlobalContext(), APInt(64, this->Displ[i], false));
        Value* inbuf_displ_int = Builder.CreateAdd(inbuf_new_int, displ_i);
        Value* inbuf_displ = Builder.CreateIntToPtr(inbuf_displ_int, Type::getInt8PtrTy(getGlobalContext()));
        Builder.CreateStore(inbuf_displ, inbuf_ptr);
        this->Types[i]->Codegen_Pack(inbuf_ptr, ConstantInt::get(getGlobalContext(), APInt(32, this->Blocklen[i], false)), outbuf_ptr);
    }

    // Emit the step value. 
    Value* StepVal = ConstantInt::get(getGlobalContext(), APInt(32, 1, false));
    Value* NextVar = Builder.CreateAdd(Variable, StepVal, "nextvar");
                            
    // check if we are finished
    Value* EndCond = Builder.CreateICmpNE(NextVar, incount, "loopcond");

    // Create the "after loop" block and insert it.
    BasicBlock *LoopEndBB = Builder.GetInsertBlock();
    BasicBlock *AfterBB = BasicBlock::Create(getGlobalContext(), "afterloop", TheFunction);

    // Insert the conditional branch into the end of LoopEndBB.
    Builder.CreateCondBr(EndCond, LoopBB, AfterBB);
                            
    // Any new code will be inserted in AfterBB.
    Builder.SetInsertPoint(AfterBB);
                            
    // Add a new entry to the PHI node for the backedge.
    Variable->addIncoming(NextVar, LoopEndBB);

    return Constant::getNullValue(Type::getInt32Ty(getGlobalContext()));

}

Value* FARC_StructDatatype::Codegen_Unpack(Value* inbuf_ptr, Value* incount, Value* outbuf_ptr) {

    // start value for the loop over incount
    Value* CounterValue = ConstantInt::get(getGlobalContext(), APInt(32, 0, false));
    // save outbuf ptr so that we can calculate outbuf = outbuf_orig + i * extend for each iter over incount
    Value* outbuf_orig = Builder.CreateLoad(outbuf_ptr);
    Value* outbuf_orig_int = Builder.CreatePtrToInt(outbuf_orig, Type::getInt64Ty(getGlobalContext()));
    Value* extend = ConstantInt::get(getGlobalContext(), APInt(32, this->getExtend(), false));
                            
    // Make the new basic block for the loop header, inserting after current block.
    Function* TheFunction = Builder.GetInsertBlock()->getParent();
    BasicBlock* PreheaderBB = Builder.GetInsertBlock();
    BasicBlock* LoopBB = BasicBlock::Create(getGlobalContext(), "loop", TheFunction);
                            
    // Insert an explicit fall through from the current block to the LoopBB.
    Builder.CreateBr(LoopBB);
                            
    // Start insertion in LoopBB.
    Builder.SetInsertPoint(LoopBB);
                            
    // Start the PHI node with an entry for Start.
    PHINode* Variable = Builder.CreatePHI(Type::getInt32Ty(getGlobalContext()), 2, "i");
    Variable->addIncoming(CounterValue, PreheaderBB);
                            
    // Emit the body of the loop:

    //outbuf = outbuf_old + incount_i * extend
    Value* incnt_i_mul_ext_32 = Builder.CreateMul(Variable, extend);
    Value* incnt_i_mul_ext_64 = Builder.CreateZExt(incnt_i_mul_ext_32, Type::getInt64Ty(getGlobalContext()));
    Value* outbuf_new_int = Builder.CreateAdd(outbuf_orig_int, incnt_i_mul_ext_64);

    for (int i=0; i<this->Count; i++) {
        // outbuf_displ = inptr_new + displ[i]
        Value* displ_i = ConstantInt::get(getGlobalContext(), APInt(64, this->Displ[i], false));
        Value* outbuf_displ_int = Builder.CreateAdd(outbuf_new_int, displ_i);
        Value* outbuf_displ = Builder.CreateIntToPtr(outbuf_displ_int, Type::getInt8PtrTy(getGlobalContext()));
        Builder.CreateStore(outbuf_displ, outbuf_ptr);
        this->Types[i]->Codegen_Unpack(inbuf_ptr, ConstantInt::get(getGlobalContext(), APInt(32, this->Blocklen[i], false)), outbuf_ptr);
    }

    // Emit the step value. 
    Value* StepVal = ConstantInt::get(getGlobalContext(), APInt(32, 1, false));
    Value* NextVar = Builder.CreateAdd(Variable, StepVal, "nextvar");
                            
    // check if we are finished
    Value* EndCond = Builder.CreateICmpNE(NextVar, incount, "loopcond");

    // Create the "after loop" block and insert it.
    BasicBlock *LoopEndBB = Builder.GetInsertBlock();
    BasicBlock *AfterBB = BasicBlock::Create(getGlobalContext(), "afterloop", TheFunction);

    // Insert the conditional branch into the end of LoopEndBB.
    Builder.CreateCondBr(EndCond, LoopBB, AfterBB);
                            
    // Any new code will be inserted in AfterBB.
    Builder.SetInsertPoint(AfterBB);
                            
    // Add a new entry to the PHI node for the backedge.
    Variable->addIncoming(NextVar, LoopEndBB);

    return Constant::getNullValue(Type::getInt32Ty(getGlobalContext()));

}

/* Class for hindexed types */
class FARC_HIndexedDatatype : public FARC_Datatype {

    int Count;
    FARC_Datatype* Basetype;
    std::vector<int> Blocklen;
    std::vector<int> Displ;

    public:
    FARC_HIndexedDatatype(int count, int* blocklen, long* displ, FARC_Datatype* basetype);
    Value *Codegen_Pack(Value* inbuf, Value* incount, Value* outbuf);
    Value *Codegen_Unpack(Value* inbuf, Value* incount, Value* outbuf);
    int getExtend();
    int getSize();

};

FARC_HIndexedDatatype::FARC_HIndexedDatatype(int count, int* blocklen, long* displ, FARC_Datatype* basetype) {

    this->Count = count;
    this->Basetype = basetype;
    for (int i=0; i<count; i++) this->Blocklen.push_back(blocklen[i]);
    for (int i=0; i<count; i++) this->Displ.push_back(displ[i]);

}

int FARC_HIndexedDatatype::getExtend() {

    if (this->Count == 0) return 0;

    int bext = this->Basetype->getExtend();

    int ub = this->Displ[0] + bext * this->Blocklen[0];
    int lb = this->Displ[0];

    for (int i=0; i<this->Count; i++) {
        int tmp_ub = this->Displ[i] + bext * this->Blocklen[i];
        int tmp_lb = this->Displ[i];
        if (tmp_ub > ub) ub = tmp_ub;
        if (tmp_lb < lb) lb = tmp_lb;
    }

    return ub - lb;

}

int FARC_HIndexedDatatype::getSize() {

    int sum = 0;
    int bsize = this->Basetype->getSize();
    for (int i=0; i<this->Count; i++) {
        sum += bsize * this->Blocklen[i];
    }

    return sum;

}

Value* FARC_HIndexedDatatype::Codegen_Pack(Value* inbuf_ptr, Value* incount, Value* outbuf_ptr) {

    // start value for the loop over incount
    Value* CounterValue = ConstantInt::get(getGlobalContext(), APInt(32, 0, false));
    // save inbuf ptr so that we can calculate inbuf = inbuf_orig + i * extend for each iter over incount
    Value* inbuf_orig = Builder.CreateLoad(inbuf_ptr);
    Value* inbuf_orig_int = Builder.CreatePtrToInt(inbuf_orig, Type::getInt64Ty(getGlobalContext()));
    Value* extend = ConstantInt::get(getGlobalContext(), APInt(32, this->getExtend(), false));
                            
    // Make the new basic block for the loop header, inserting after current block.
    Function* TheFunction = Builder.GetInsertBlock()->getParent();
    BasicBlock* PreheaderBB = Builder.GetInsertBlock();
    BasicBlock* LoopBB = BasicBlock::Create(getGlobalContext(), "loop", TheFunction);
                            
    // Insert an explicit fall through from the current block to the LoopBB.
    Builder.CreateBr(LoopBB);
                            
    // Start insertion in LoopBB.
    Builder.SetInsertPoint(LoopBB);
                            
    // Start the PHI node with an entry for Start.
    PHINode* Variable = Builder.CreatePHI(Type::getInt32Ty(getGlobalContext()), 2, "i");
    Variable->addIncoming(CounterValue, PreheaderBB);
                            
    // Emit the body of the loop:

    //inbuf = inbuf_old + incount_i * extend
    Value* incnt_i_mul_ext_32 = Builder.CreateMul(Variable, extend);
    Value* incnt_i_mul_ext_64 = Builder.CreateZExt(incnt_i_mul_ext_32, Type::getInt64Ty(getGlobalContext()));
    Value* inbuf_new_int = Builder.CreateAdd(inbuf_orig_int, incnt_i_mul_ext_64);

    for (int i=0; i<this->Count; i++) {
        // inbuf_displ = inptr_new + displ[i]
        Value* displ_i = ConstantInt::get(getGlobalContext(), APInt(64, this->Displ[i], false));
        Value* inbuf_displ_int = Builder.CreateAdd(inbuf_new_int, displ_i);
        Value* inbuf_displ = Builder.CreateIntToPtr(inbuf_displ_int, Type::getInt8PtrTy(getGlobalContext()));
        Builder.CreateStore(inbuf_displ, inbuf_ptr);
        Basetype->Codegen_Pack(inbuf_ptr, ConstantInt::get(getGlobalContext(), APInt(32, this->Blocklen[i], false)), outbuf_ptr);
    }

    // Emit the step value. 
    Value* StepVal = ConstantInt::get(getGlobalContext(), APInt(32, 1, false));
    Value* NextVar = Builder.CreateAdd(Variable, StepVal, "nextvar");
                            
    // check if we are finished
    Value* EndCond = Builder.CreateICmpNE(NextVar, incount, "loopcond");

    // Create the "after loop" block and insert it.
    BasicBlock *LoopEndBB = Builder.GetInsertBlock();
    BasicBlock *AfterBB = BasicBlock::Create(getGlobalContext(), "afterloop", TheFunction);

    // Insert the conditional branch into the end of LoopEndBB.
    Builder.CreateCondBr(EndCond, LoopBB, AfterBB);
                            
    // Any new code will be inserted in AfterBB.
    Builder.SetInsertPoint(AfterBB);
                            
    // Add a new entry to the PHI node for the backedge.
    Variable->addIncoming(NextVar, LoopEndBB);

    return Constant::getNullValue(Type::getInt32Ty(getGlobalContext()));


}

Value* FARC_HIndexedDatatype::Codegen_Unpack(Value* inbuf_ptr, Value* incount, Value* outbuf_ptr) {

    // start value for the loop over incount
    Value* CounterValue = ConstantInt::get(getGlobalContext(), APInt(32, 0, false));
    // save inbuf ptr so that we can calculate inbuf = inbuf_orig + i * extend for each iter over incount
    Value* outbuf_orig = Builder.CreateLoad(outbuf_ptr);
    Value* outbuf_orig_int = Builder.CreatePtrToInt(outbuf_orig, Type::getInt64Ty(getGlobalContext()));
    Value* extend = ConstantInt::get(getGlobalContext(), APInt(32, this->getExtend(), false));
                            
    // Make the new basic block for the loop header, inserting after current block.
    Function* TheFunction = Builder.GetInsertBlock()->getParent();
    BasicBlock* PreheaderBB = Builder.GetInsertBlock();
    BasicBlock* LoopBB = BasicBlock::Create(getGlobalContext(), "loop", TheFunction);
                            
    // Insert an explicit fall through from the current block to the LoopBB.
    Builder.CreateBr(LoopBB);
                            
    // Start insertion in LoopBB.
    Builder.SetInsertPoint(LoopBB);
                            
    // Start the PHI node with an entry for Start.
    PHINode* Variable = Builder.CreatePHI(Type::getInt32Ty(getGlobalContext()), 2, "i");
    Variable->addIncoming(CounterValue, PreheaderBB);
                            
    // Emit the body of the loop:

    //inbuf = inbuf_old + incount_i * extend
    Value* incnt_i_mul_ext_32 = Builder.CreateMul(Variable, extend);
    Value* incnt_i_mul_ext_64 = Builder.CreateZExt(incnt_i_mul_ext_32, Type::getInt64Ty(getGlobalContext()));
    Value* outbuf_new_int = Builder.CreateAdd(outbuf_orig_int, incnt_i_mul_ext_64);

    for (int i=0; i<this->Count; i++) {
        // inbuf_displ = inptr_new + displ[i]
        Value* displ_i = ConstantInt::get(getGlobalContext(), APInt(64, this->Displ[i], false));
        Value* outbuf_displ_int = Builder.CreateAdd(outbuf_new_int, displ_i);
        Value* outbuf_displ = Builder.CreateIntToPtr(outbuf_displ_int, Type::getInt8PtrTy(getGlobalContext()));
        Builder.CreateStore(outbuf_displ, outbuf_ptr);
        Basetype->Codegen_Unpack(inbuf_ptr, ConstantInt::get(getGlobalContext(), APInt(32, this->Blocklen[i], false)), outbuf_ptr);
    }

    // Emit the step value. 
    Value* StepVal = ConstantInt::get(getGlobalContext(), APInt(32, 1, false));
    Value* NextVar = Builder.CreateAdd(Variable, StepVal, "nextvar");
                            
    // check if we are finished
    Value* EndCond = Builder.CreateICmpNE(NextVar, incount, "loopcond");

    // Create the "after loop" block and insert it.
    BasicBlock *LoopEndBB = Builder.GetInsertBlock();
    BasicBlock *AfterBB = BasicBlock::Create(getGlobalContext(), "afterloop", TheFunction);

    // Insert the conditional branch into the end of LoopEndBB.
    Builder.CreateCondBr(EndCond, LoopBB, AfterBB);
                            
    // Any new code will be inserted in AfterBB.
    Builder.SetInsertPoint(AfterBB);
                            
    // Add a new entry to the PHI node for the backedge.
    Variable->addIncoming(NextVar, LoopEndBB);

    return Constant::getNullValue(Type::getInt32Ty(getGlobalContext()));


}

// this jits the pack/unpack functions
void generate_pack_function(FARC_Datatype* ddt) {

    std::vector<std::string> Args;
    Args.push_back("inbuf");
    Args.push_back("count");
    Args.push_back("outbuf");

    std::vector<Type*> FuncArgs;
    FuncArgs.push_back(Type::getInt8PtrTy(getGlobalContext()));
    FuncArgs.push_back(Type::getInt32Ty(getGlobalContext()));
    FuncArgs.push_back(Type::getInt8PtrTy(getGlobalContext()));

    FunctionType *FT = FunctionType::get(Type::getInt32Ty(getGlobalContext()), FuncArgs, false);

    Function* F = Function::Create(FT, Function::ExternalLinkage, "packer", TheModule);
    // Set names for all arguments.
    unsigned Idx = 0;

    for (Function::arg_iterator AI = F->arg_begin(); Idx != Args.size(); ++AI, ++Idx) {
        AI->setName(Args[Idx]);
        NamedValues[Args[Idx]] = AI;
    }

    // Create a new basic block to start insertion into.
    BasicBlock *BB = BasicBlock::Create(getGlobalContext(), "entry", F);
    Builder.SetInsertPoint(BB);

    // generate code for the datatype
    Value* RetVal = ddt->Codegen_Pack(NamedValues["inbuf"], NamedValues["count"], NamedValues["outbuf"]);
    Builder.CreateRet(RetVal);

//    verifyFunction(*F);
//    TheFPM->run(*F);
//    F->dump();

    TheExecutionEngine->runJITOnFunction(F);
    ddt->packer = TheExecutionEngine->getPointerToFunction(F);

}

void generate_unpack_function(FARC_Datatype* ddt) {

    std::vector<std::string> Args;
    Args.push_back("inbuf");
    Args.push_back("count");
    Args.push_back("outbuf");

    std::vector<Type*> FuncArgs;
    FuncArgs.push_back(Type::getInt8PtrTy(getGlobalContext()));
    FuncArgs.push_back(Type::getInt32Ty(getGlobalContext()));
    FuncArgs.push_back(Type::getInt8PtrTy(getGlobalContext()));

    FunctionType *FT = FunctionType::get(Type::getInt32Ty(getGlobalContext()), FuncArgs, false);

    Function* F = Function::Create(FT, Function::ExternalLinkage, "unpacker", TheModule);
    // Set names for all arguments.
    unsigned Idx = 0;

    for (Function::arg_iterator AI = F->arg_begin(); Idx != Args.size(); ++AI, ++Idx) {
        AI->setName(Args[Idx]);
        NamedValues[Args[Idx]] = AI;
    }

    // Create a new basic block to start insertion into.
    BasicBlock *BB = BasicBlock::Create(getGlobalContext(), "entry", F);
    Builder.SetInsertPoint(BB);

    // generate code for the datatype
    Value* RetVal = ddt->Codegen_Unpack(NamedValues["inbuf"], NamedValues["count"], NamedValues["outbuf"]);
    Builder.CreateRet(RetVal);

//    verifyFunction(*F);
//    TheFPM->run(*F);
    F->dump();

    TheExecutionEngine->runJITOnFunction(F);
    ddt->unpacker = TheExecutionEngine->getPointerToFunction(F);

}

FARC_Datatype* FARC_DDT_Commit(FARC_Datatype* ddt) {

//    This is usefull for debugging, as it allows us to call printf during jitting
//    std::vector<llvm::Type*> printf_arg_types;
//    printf_arg_types.push_back(Type::getInt8PtrTy(getGlobalContext()));
//    FunctionType* printf_type = FunctionType::get(Type::getInt32Ty(getGlobalContext()), printf_arg_types, true);
//    func = Function::Create(printf_type, Function::ExternalLinkage, Twine("printf"), TheModule);

    generate_pack_function(ddt);
    generate_unpack_function(ddt);

    return ddt;

}

// this calls the pack/unpack function
void FARC_DDT_Pack(char* inbuf, char* outbuf, FARC_Datatype* ddt, int count) {

     int (*FP)(void*, int, void*) = (int (*)(void*, int, void*))(intptr_t) ddt->packer;
     FP(inbuf, count, outbuf);
    
}

void FARC_DDT_Unpack(char* inbuf, char* outbuf, FARC_Datatype* ddt, int count) {

     int (*FP)(void*, int, void*) = (int (*)(void*, int, void*))(intptr_t) ddt->unpacker;
     FP(inbuf, count, outbuf);
    
}

void FARC_DDT_Free(FARC_Datatype* ddt) {

    delete ddt;
    
}

// init the JIT compiler
void FARC_DDT_Init() {

    InitializeNativeTarget();
    LLVMContext &Context = getGlobalContext();
    TheModule = new Module("FARC-JIT", Context);


    // Create the JIT.  This takes ownership of the module.
    std::string ErrStr;
    TheExecutionEngine = EngineBuilder(TheModule).setErrorStr(&ErrStr).create();
    if (!TheExecutionEngine) {
        fprintf(stderr, "Could not create ExecutionEngine: %s\n", ErrStr.c_str());
        exit(1);
    }

/*
    FunctionPassManager* OurFPM = new FunctionPassManager(TheModule);

//    PassManagerBuilder Builder;
//    Builder.OptLevel = 3;
//    Builder.Vectorize = true;
//    Builder.LoopVectorize = true;
//    Builder.populateFunctionPassManager(*OurFPM);


    // Set up the optimizer pipeline.  Start with registering info about how the
    // target lays out data structures.
    OurFPM->add(new DataLayout(*TheExecutionEngine->getDataLayout()));

    OurFPM->add(createBasicAliasAnalysisPass());  // -basicaa
    OurFPM->add(createPromoteMemoryToRegisterPass()); // -mem2reg
    OurFPM->add(createCFGSimplificationPass());   // -simplifycfg
    OurFPM->add(createInstructionCombiningPass());    // -instcombine
//    OurFPM->add(createReassociatePass());
    OurFPM->add(createGVNPass());
    OurFPM->add(createCFGSimplificationPass());
//    OurFPM->add(createTailCallEliminationPass()); // -tailcallelim
//    OurFPM->add(createLoopSimplifyPass());        // -loop-simplify
//    OurFPM->add(createLCSSAPass());           // -lcssa
//    OurFPM->add(createLoopRotatePass());      // -loop-rotate
    OurFPM->add(createLCSSAPass());           // -lcssa
//    OurFPM->add(createLoopUnswitchPass());        // -loop-unswitch
//    OurFPM->add(createInstructionCombiningPass());    // -instcombine
    OurFPM->add(createLoopSimplifyPass());        // -loop-simplify
//    OurFPM->add(createLCSSAPass());           // -lcssa
    OurFPM->add(createIndVarSimplifyPass());      // -indvars
    OurFPM->add(createLoopUnrollPass());
//    OurFPM->add(createLoopDeletionPass());        // -loop-deletion
    OurFPM->add(createInstructionCombiningPass());    // -instcombine
//    OurFPM->add(createLoopVectorizePass());
//    OurFPM->add(createBBVectorizePass());

    OurFPM->add(createAggressiveDCEPass());
    OurFPM->doInitialization();


    // Set the global so the code gen can use this.
    TheFPM = OurFPM;
*/

}

#endif

