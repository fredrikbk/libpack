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


/* Codegen Functions */
Value* multNode(int op1, Value* op2PtrNode) {
    Value* op1Node = ConstantInt::get(getGlobalContext(), APInt(64, op1, false));
    Value* op2Node = Builder.CreateIntCast(op2PtrNode, Type::getInt64Ty(getGlobalContext()), false); 
    return Builder.CreateMul(op1Node, op2Node);
}

Value* constNode(int val) {
    return ConstantInt::get(getGlobalContext(), APInt(32, val, false));
}

Value* constNode(long val) {
    return ConstantInt::get(getGlobalContext(), APInt(64, val, false));
}

void vectorCodegen(Value* inbuf, Value* incount, Value* outbuf, FARC_Datatype* basetype,
        int count, int blocklen, int elemstride_in, int elemstride_out, bool pack) {
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
    basetype->Codegen_Pack(in_inner, ConstantInt::get(getGlobalContext(), APInt(32, blocklen, false)), out_inner);


    // Increment the out ptr by Size(Basetype) * Blocklen
    Value* out_bytes_to_stride = ConstantInt::get(getGlobalContext(), APInt(64, elemstride_out, false));
    Value* out_addr_cvi = Builder.CreatePtrToInt(out_inner, Type::getInt64Ty(getGlobalContext()));
    Value* out_addr = Builder.CreateAdd(out_addr_cvi, out_bytes_to_stride);
    Value* nextout_inner = Builder.CreateIntToPtr(out_addr, Type::getInt8PtrTy(getGlobalContext()));

    // Increment the in ptr by Extend(Basetype) * Stride
    Value* in_bytes_to_stride = ConstantInt::get(getGlobalContext(), APInt(64, elemstride_in, false));
    Value* in_addr_cvi = Builder.CreatePtrToInt(in_inner, Type::getInt64Ty(getGlobalContext()));
    Value* in_addr = Builder.CreateAdd(in_addr_cvi, in_bytes_to_stride);
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
    Value* EndCond_inner = Builder.CreateICmpEQ(nextj, constNode(count), "innercond");
                            
    // Create and branch to the inner loop postamble
    BasicBlock *LoopEnd_inner_BB = Builder.GetInsertBlock();
    BasicBlock *After_inner_BB = BasicBlock::Create(getGlobalContext(), "afterinner", TheFunction);
    Builder.CreateCondBr(EndCond_inner, After_inner_BB, Loop_inner_BB);
    Builder.SetInsertPoint(After_inner_BB);

    // Add backedges for the inner loop induction variables
    out_inner->addIncoming(nextout_inner, LoopEnd_inner_BB);
    in_inner->addIncoming(nextin_inner, LoopEnd_inner_BB);
    j->addIncoming(nextj, LoopEnd_inner_BB);


    // Move the the extend-stride ptr back Extend(Basetype) * Stride - Size(Basetype) * Blocklen  
    Value* nextin_outer = NULL;
    Value* nextout_outer = NULL;
    if (pack) {
        nextout_outer = nextout_inner;
        Value* nextin_outer_val = Builder.CreateAdd(in_addr_cvi, out_bytes_to_stride );
        nextin_outer = Builder.CreateIntToPtr(nextin_outer_val, Type::getInt8PtrTy(getGlobalContext()));
    }
    else {
        Value* nextout_outer_val = Builder.CreateAdd(out_addr_cvi, in_bytes_to_stride );
        nextout_outer = Builder.CreateIntToPtr(nextout_outer_val, Type::getInt8PtrTy(getGlobalContext()));
        nextin_outer = nextin_inner; 
    }

    // Increment outer loop index
    Value* nexti = Builder.CreateAdd(i, constNode(1), "nexti");
    Value* EndCond_outer = Builder.CreateICmpEQ(nexti, incount, "outercond");

    // Create and branch to the outer loop postamble
    BasicBlock *LoopEnd_outer_BB = Builder.GetInsertBlock();
    BasicBlock *After_outer_BB = BasicBlock::Create(getGlobalContext(), "afterouter", TheFunction);
    Builder.CreateCondBr(EndCond_outer, After_outer_BB, Loop_outer_BB);
    Builder.SetInsertPoint(After_outer_BB);

    // Add backedges for the outer loop induction variable
    out_outer->addIncoming(nextout_outer, LoopEnd_outer_BB);
    in_outer->addIncoming(nextin_outer, LoopEnd_outer_BB);
    i->addIncoming(nexti, LoopEnd_outer_BB);
}


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

Value* FARC_PrimitiveDatatype::Codegen_Unpack(Value* inbuf, Value* incount, Value* outbuf) {

    // does exactly the same as pack for primitive types
    return this->Codegen_Pack(inbuf, incount, outbuf);
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
    void Codegen(Value* inbuf, Value* incount, Value* outbuf, int elemstride_in, int elemstride_out);

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

void FARC_ContiguousDatatype::Codegen(Value* inbuf, Value* incount, Value* outbuf, int elemstride_in, int elemstride_out) {
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
    Value* out_bytes_to_stride = ConstantInt::get(getGlobalContext(), APInt(64, elemstride_out * this->Count, false));
    Value* out_addr_cvi = Builder.CreatePtrToInt(out, Type::getInt64Ty(getGlobalContext()));
    Value* out_addr = Builder.CreateAdd(out_addr_cvi, out_bytes_to_stride);
    Value* nextout = Builder.CreateIntToPtr(out_addr, Type::getInt8PtrTy(getGlobalContext()));

    // Increment the in ptr by Extend(Basetype) * Stride
    Value* in_bytes_to_stride = ConstantInt::get(getGlobalContext(), APInt(64, elemstride_in * this->Count, false));
    Value* in_addr_cvi = Builder.CreatePtrToInt(in, Type::getInt64Ty(getGlobalContext()));
    Value* in_addr = Builder.CreateAdd(in_addr_cvi, in_bytes_to_stride);
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
}

Value* FARC_ContiguousDatatype::Codegen_Pack(Value* inbuf, Value* incount, Value* outbuf) {
    Codegen(inbuf, incount, outbuf, this->Basetype->getExtend(), this->Basetype->getSize());
    return Constant::getNullValue(Type::getInt32Ty(getGlobalContext()));
}

Value* FARC_ContiguousDatatype::Codegen_Unpack(Value* inbuf, Value* incount, Value* outbuf) {
    Codegen(inbuf, incount, outbuf, this->Basetype->getSize(), this->Basetype->getExtend());
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

Value* FARC_HVectorDatatype::Codegen_Pack(Value* inbuf, Value* incount, Value* outbuf) {
    vectorCodegen(inbuf, incount, outbuf, this->Basetype, this->Count,
            this->Blocklen, this->Stride, this->Basetype->getSize() * this->Blocklen, true);
    return Constant::getNullValue(Type::getInt32Ty(getGlobalContext()));
}

Value* FARC_HVectorDatatype::Codegen_Unpack(Value* inbuf, Value* incount, Value* outbuf) {
    vectorCodegen(inbuf, incount, outbuf, this->Basetype, this->Count, this->Blocklen,
            this->Basetype->getSize() * this->Blocklen, this->Stride, false);
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
    vectorCodegen(inbuf, incount, outbuf, this->Basetype, this->Count,
            this->Blocklen, this->Basetype->getExtend() * this->Stride, this->Basetype->getSize() * this->Blocklen, true);
    return Constant::getNullValue(Type::getInt32Ty(getGlobalContext()));
}

Value* FARC_VectorDatatype::Codegen_Unpack(Value* inbuf, Value* incount, Value* outbuf) {
    vectorCodegen(inbuf, incount, outbuf, this->Basetype, this->Count, this->Blocklen,
            this->Basetype->getSize() * this->Blocklen, this->Basetype->getExtend() * this->Stride, false);
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
    void Codegen(Value *contig_buf, Value *hindexed_buf, Value* incount, bool gather);

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

void FARC_HIndexedDatatype::Codegen(Value *compactbuf, Value *scatteredbuf, Value* incount, bool pack) {
    Function* TheFunction = Builder.GetInsertBlock()->getParent();

    // Base address of the input buffer
    Value* scatteredbuf_orig_int = Builder.CreatePtrToInt(scatteredbuf, Type::getInt64Ty(getGlobalContext()));
    Value* extend = constNode((long)this->getExtend());
    Value* incount_64 = Builder.CreateZExt(incount, Type::getInt64Ty(getGlobalContext()));
    Value* incount_expanded = Builder.CreateMul(incount_64, extend);

    // Loop
    BasicBlock* PreheaderBB = Builder.GetInsertBlock();
    BasicBlock* LoopBB = BasicBlock::Create(getGlobalContext(), "loop", TheFunction);
    Builder.CreateBr(LoopBB);
    Builder.SetInsertPoint(LoopBB);

    PHINode *compact = Builder.CreatePHI(Type::getInt8PtrTy(getGlobalContext()), 2, "compact");
    compact->addIncoming(compactbuf, PreheaderBB);
    PHINode* i = Builder.CreatePHI(Type::getInt64Ty(getGlobalContext()), 2, "i");
    i->addIncoming(constNode(0l), PreheaderBB);

    Value* compact_addr = Builder.CreatePtrToInt(compact, Type::getInt64Ty(getGlobalContext()));

    // OPT: Make this the loop counter
    Value* scattered_disp_base = Builder.CreateAdd(scatteredbuf_orig_int, i);

    Value* nextcompact = compact;
    for (int i=0; i<this->Count; i++) {
        // Set the scattered ptr to scattered_disp_base + this->Disl[i]
        Value* displ_i = ConstantInt::get(getGlobalContext(), APInt(64, this->Displ[i], false));
        Value* scattered_disp = Builder.CreateAdd(scattered_disp_base, displ_i);
        Value* scattered = Builder.CreateIntToPtr(scattered_disp, Type::getInt8PtrTy(getGlobalContext()));

        if (pack) {
            Basetype->Codegen_Pack(scattered, ConstantInt::get(getGlobalContext(), APInt(32, this->Blocklen[i], false)), nextcompact);
        }
        else {
            Basetype->Codegen_Unpack(nextcompact, ConstantInt::get(getGlobalContext(), APInt(32, this->Blocklen[i], false)), scattered);
        }

        // Increment the compact ptr by Size(Basetype) * Blocklen
        Value* compact_bytes_to_stride = constNode((long)Basetype->getSize() * this->Blocklen[i]);
        compact_addr = Builder.CreateAdd(compact_addr, compact_bytes_to_stride);
        nextcompact = Builder.CreateIntToPtr(compact_addr, Type::getInt8PtrTy(getGlobalContext()));
    }

    // Increment the loop index and test for loop exit
    Value* nexti = Builder.CreateAdd(i, extend, "nexti");
    Value* EndCond = Builder.CreateICmpEQ(nexti, incount_expanded, "loopcond");

    // Create and branch to the outer loop postamble
    BasicBlock *LoopEndBB = Builder.GetInsertBlock();
    BasicBlock *AfterBB = BasicBlock::Create(getGlobalContext(), "afterloop", TheFunction);

    Builder.CreateCondBr(EndCond, AfterBB, LoopBB);
    Builder.SetInsertPoint(AfterBB);
                            
    // Add backedges for the loop induction variable
    compact->addIncoming(nextcompact, LoopEndBB);
    i->addIncoming(nexti, LoopEndBB);
}

Value* FARC_HIndexedDatatype::Codegen_Pack(Value* inbuf, Value* incount, Value* outbuf) {
    Codegen(outbuf, inbuf, incount, true);
    return Constant::getNullValue(Type::getInt32Ty(getGlobalContext()));
}

Value* FARC_HIndexedDatatype::Codegen_Unpack(Value* inbuf, Value* incount, Value* outbuf) {
    Codegen(inbuf, outbuf, incount, false);
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
//    F->dump();

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

