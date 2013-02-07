#ifndef DDT_JIT_HPP
#define DDT_JIT_HPP

#include <map>
#include <string>
#include <vector>
#include <cstdio>

#include <mpi.h>

#include "llvm/IR/Module.h"
#include "llvm/IR/LLVMContext.h"
#include "llvm/IR/DerivedTypes.h"
#include "llvm/IR/IRBuilder.h"
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

static std::vector<void*> Packers;
static Module *TheModule;
static IRBuilder<> Builder(getGlobalContext());
static std::map<std::string, Value*> NamedValues;
static ExecutionEngine *TheExecutionEngine;
static FunctionPassManager *TheFPM;



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

Value* FARC_PrimitiveDatatype::Codegen_Pack(Value* inbuf_ptr, Value* incount, Value* outbuf_ptr) {
                            
    /** for debugging  
    std::vector<llvm::Type*> printf_arg_types;
    printf_arg_types.push_back(Type::getInt8PtrTy(getGlobalContext()));
    FunctionType* printf_type = FunctionType::get(Type::getInt32Ty(getGlobalContext()), printf_arg_types, true);
    Function *func = Function::Create(printf_type, Function::ExternalLinkage, Twine("printf"), TheModule);
    Value *fmt_ptr1 = Builder.CreateGlobalStringPtr("Before in PTR: %p Out PTR: %p\n\0");
    // now we can print as follows:
    CallInst *call = Builder.CreateCall3(func, fmt_ptr1, inbuf, outbuf);
*/

    Value* inbuf = Builder.CreateLoad(inbuf_ptr);
    Value* outbuf = Builder.CreateLoad(outbuf_ptr);

    Value* contig_size_type = ConstantInt::get(getGlobalContext(), APInt(64, this->Extend, false));
    Value* incount64 = Builder.CreateIntCast(incount, Type::getInt64Ty(getGlobalContext()), false); 
    Value* contig_size = Builder.CreateMul(contig_size_type, incount64);

    Value* memcopy = Builder.CreateMemCpy(outbuf, inbuf, contig_size, 1);
//    CallInst *call = Builder.CreateCall4(func, fmt_ptr1, inbuf, outbuf, contig_size);

    Value* outbuf_addr_cvi = Builder.CreatePtrToInt(outbuf, Type::getInt64Ty(getGlobalContext()));
    Value* outbuf_addr = Builder.CreateAdd(outbuf_addr_cvi, contig_size);
    Value* outbuf_addr_cvp = Builder.CreateIntToPtr(outbuf_addr, Type::getInt8PtrTy(getGlobalContext()));
                            
    Value* inbuf_addr_cvi = Builder.CreatePtrToInt(inbuf, Type::getInt64Ty(getGlobalContext()));
    Value* inbuf_addr = Builder.CreateAdd(inbuf_addr_cvi, contig_size);
    Value* inbuf_addr_cvp = Builder.CreateIntToPtr(inbuf_addr, Type::getInt8PtrTy(getGlobalContext()));

    Builder.CreateStore(inbuf_addr_cvp, inbuf_ptr);
    Builder.CreateStore(outbuf_addr_cvp, outbuf_ptr);

#if 0                            
    Value* inbuf = Builder.CreateLoad(inbuf_ptr);
    Value* outbuf = Builder.CreateLoad(outbuf_ptr);

    // start value for the loop over incount
    Value* CounterValue = ConstantInt::get(getGlobalContext(), APInt(64, 0, false));
                            
    // end value            
    Value* MPI_Type_size = ConstantInt::get(getGlobalContext(), APInt(32, this->Extend, false));
    Value* EndVal32 = Builder.CreateMul(incount, MPI_Type_size);
    Value* EndVal = Builder.CreateZExt(EndVal32, Type::getInt64Ty(getGlobalContext()));
                            
    // Make the new basic block for the loop header, inserting after current block.
    Function* TheFunction = Builder.GetInsertBlock()->getParent();
    BasicBlock* PreheaderBB = Builder.GetInsertBlock();
    BasicBlock* LoopBB = BasicBlock::Create(getGlobalContext(), "loop", TheFunction);
                            
    // Insert an explicit fall through from the current block to the LoopBB.
    Builder.CreateBr(LoopBB);
                            
    // Start insertion in LoopBB.
    Builder.SetInsertPoint(LoopBB);
                            
    // Start the PHI node with an entry for Start.
    PHINode* Variable = Builder.CreatePHI(Type::getInt64Ty(getGlobalContext()), 2, "i");
    Variable->addIncoming(CounterValue, PreheaderBB);
                            
    // Emit the body of the loop.
                            
 //    CallInst *call = Builder.CreateCall3(func, fmt_ptr, inbuf_addr_cvp, outbuf_addr_cvp);
                            
   Value* outbuf_addr_cvi = Builder.CreatePtrToInt(outbuf, Type::getInt64Ty(getGlobalContext()));
   Value* outbuf_addr = Builder.CreateAdd(outbuf_addr_cvi, Variable);
   Value* outbuf_addr_cvp = Builder.CreateIntToPtr(outbuf_addr, Type::getInt8PtrTy(getGlobalContext()));
                            
   Value* inbuf_addr_cvi = Builder.CreatePtrToInt(inbuf, Type::getInt64Ty(getGlobalContext()));
   Value* inbuf_addr = Builder.CreateAdd(inbuf_addr_cvi, Variable);
   Value* inbuf_addr_cvp = Builder.CreateIntToPtr(inbuf_addr, Type::getInt8PtrTy(getGlobalContext()));

    Value* tmp = Builder.CreateLoad(inbuf_addr_cvp);
    Builder.CreateStore(tmp, outbuf_addr_cvp);
                            
    Builder.CreateStore(inbuf_addr_cvp, inbuf_ptr);
    Builder.CreateStore(outbuf_addr_cvp, outbuf_ptr);

    // Emit the step value. 
    Value* StepVal = ConstantInt::get(getGlobalContext(), APInt(64, 1, false));
    Value* NextVar = Builder.CreateAdd(Variable, StepVal, "nextvar");
                            
    // check if we are finished
    Value* EndCond = Builder.CreateICmpNE(NextVar, EndVal, "loopcond");

    // Create the "after loop" block and insert it.
    BasicBlock *LoopEndBB = Builder.GetInsertBlock();
    BasicBlock *AfterBB = BasicBlock::Create(getGlobalContext(), "afterloop", TheFunction);

    // Insert the conditional branch into the end of LoopEndBB.
    Builder.CreateCondBr(EndCond, LoopBB, AfterBB);
                            
    // Any new code will be inserted in AfterBB.
    Builder.SetInsertPoint(AfterBB);
                            
   Value* outbuf_addr2 = Builder.CreateAdd(outbuf_addr_cvi, EndVal);
   Value* outbuf_addr_cvp2 = Builder.CreateIntToPtr(outbuf_addr2, Type::getInt8PtrTy(getGlobalContext()));
                            
   Value* inbuf_addr2 = Builder.CreateAdd(inbuf_addr_cvi, EndVal);
   Value* inbuf_addr_cvp2 = Builder.CreateIntToPtr(inbuf_addr2, Type::getInt8PtrTy(getGlobalContext()));
                            
   Builder.CreateStore(inbuf_addr_cvp2, inbuf_ptr);
   Builder.CreateStore(outbuf_addr_cvp2, outbuf_ptr);

    // Add a new entry to the PHI node for the backedge.
    Variable->addIncoming(NextVar, LoopEndBB);
#endif

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

Value* FARC_ContiguousDatatype::Codegen_Pack(Value* inbuf_ptr, Value* incount, Value* outbuf_ptr) {

    // start value for the loop over incount
    Value* CounterValue = ConstantInt::get(getGlobalContext(), APInt(32, 0, false));
                            
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
                            
    // Emit the body of the loop.
                            
    // Codegen for subtype                            
    Basetype->Codegen_Pack(inbuf_ptr, ConstantInt::get(getGlobalContext(), APInt(32, this->Count, false)), outbuf_ptr);
                            
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

Value* FARC_ContiguousDatatype::Codegen_Unpack(Value* inbuf_ptr, Value* incount, Value* outbuf_ptr) {

    return this->Codegen_Pack(inbuf_ptr, incount, outbuf_ptr);

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
                            
    Basetype->Codegen_Pack(inbuf_ptr, ConstantInt::get(getGlobalContext(), APInt(32, this->Blocklen, false)), outbuf_ptr);
                            
    Value* newout = Builder.CreateLoad(outbuf_ptr);
                            
    // *imbuf += stride - (newin - oldin)
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
    Builder.CreateStore(newout, inbuf_ptr);
                            
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

Value* FARC_VectorDatatype::Codegen_Pack(Value* inbuf_ptr, Value* incount, Value* outbuf_ptr) {

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
    Value* stride = ConstantInt::get(getGlobalContext(), APInt(64, Stride*this->Basetype->getExtend(), false));
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
                            
    Basetype->Codegen_Pack(inbuf_ptr, ConstantInt::get(getGlobalContext(), APInt(32, this->Blocklen, false)), outbuf_ptr);
                            
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
    Builder.CreateStore(newout, inbuf_ptr);
                            
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

// this jits the pack/unpack functions
int FARC_DDT_Commit(FARC_Datatype* ddt) {

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

    Value* inbuf_ptr = Builder.CreateAlloca(Type::getInt8PtrTy(getGlobalContext()));
    Value* outbuf_ptr = Builder.CreateAlloca(Type::getInt8PtrTy(getGlobalContext()));

    Builder.CreateStore(NamedValues["inbuf"], inbuf_ptr);
    Builder.CreateStore(NamedValues["outbuf"], outbuf_ptr);

    // generate code for the datatype
    Value* RetVal = ddt->Codegen_Pack(inbuf_ptr, NamedValues["count"], outbuf_ptr);

    Builder.CreateRet(RetVal);

    verifyFunction(*F);
    TheFPM->run(*F);

    //F->dump();

     // store a function pointer to the JITed function in a lookup table,
     // return lookup table index
     void *FPtr = TheExecutionEngine->getPointerToFunction(F);

     Packers.push_back(FPtr);

     ddt->packer = FPtr;

     return Packers.size()-1;

}

// this calls the pack/unpack function
void FARC_DDT_Pack(char* inbuf, char* outbuf, int ddt_handle, int count) {

     int (*FP)(void*, int, void*) = (int (*)(void*, int, void*))(intptr_t) Packers[ddt_handle];
     FP(inbuf, count, outbuf);
    
}


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

    FunctionPassManager* OurFPM = new FunctionPassManager(TheModule);

    PassManagerBuilder Builder;
    Builder.OptLevel = 3;
    Builder.Vectorize = true;
    Builder.LoopVectorize = true;
    Builder.populateFunctionPassManager(*OurFPM);

/*
    // Set up the optimizer pipeline.  Start with registering info about how the
    // target lays out data structures.
    OurFPM->add(new DataLayout(*TheExecutionEngine->getDataLayout()));

    OurFPM->add(createBasicAliasAnalysisPass());  // -basicaa
    OurFPM->add(createPromoteMemoryToRegisterPass()); // -mem2reg
    OurFPM->add(createCFGSimplificationPass());   // -simplifycfg
    OurFPM->add(createInstructionCombiningPass());    // -instcombine
    OurFPM->add(createReassociatePass());
    OurFPM->add(createGVNPass());
    OurFPM->add(createCFGSimplificationPass());
    OurFPM->add(createTailCallEliminationPass()); // -tailcallelim
    OurFPM->add(createLoopSimplifyPass());        // -loop-simplify
    OurFPM->add(createLCSSAPass());           // -lcssa
    OurFPM->add(createLoopRotatePass());      // -loop-rotate
    OurFPM->add(createLCSSAPass());           // -lcssa
    OurFPM->add(createLoopUnswitchPass());        // -loop-unswitch
    OurFPM->add(createInstructionCombiningPass());    // -instcombine
    OurFPM->add(createLoopSimplifyPass());        // -loop-simplify
    OurFPM->add(createLCSSAPass());           // -lcssa
    OurFPM->add(createIndVarSimplifyPass());      // -indvars
    OurFPM->add(createLoopDeletionPass());        // -loop-deletion
    OurFPM->add(createInstructionCombiningPass());    // -instcombine
*/

    OurFPM->doInitialization();


    // Set the global so the code gen can use this.
    TheFPM = OurFPM;

}

#endif

