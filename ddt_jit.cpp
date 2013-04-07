
#include "ddt_jit.hpp"
#include "pack.h"

#include <map>
#include <cstdio>
#include <sstream>

#if defined(LLVM_VERSION_MINOR)
#error fpp
#endif

#if defined(LLVM_VERSION_MINOR) && LLVM_VERSION_MINOR < 3
#define LLVM32 1
#else
#define LLVM32 0
#endif


#if LLVM32
#include "llvm/Module.h"
#include "llvm/LLVMContext.h"
#include "llvm/DerivedTypes.h"
#include "llvm/IRBuilder.h"
#include "llvm/Intrinsics.h"
#else
#include "llvm/IR/Module.h"
#include "llvm/IR/LLVMContext.h"
#include "llvm/IR/DerivedTypes.h"
#include "llvm/IR/IRBuilder.h"
#include "llvm/IR/Intrinsics.h"
#endif

#include "llvm/Support/TargetSelect.h"
#include "llvm/ExecutionEngine/ExecutionEngine.h"
#include "llvm/ExecutionEngine/JIT.h"

#define LLVM_VERIFY LLVM_OUTPUT
#if LLVM_VERIFY
#include "llvm/Analysis/Verifier.h"
#endif

#if LLVM_OPTIMIZE
#include "llvm/PassManager.h"
#include "llvm/Analysis/Passes.h"
#include "llvm/LinkAllPasses.h"
#include "llvm/Transforms/Scalar.h"
#include "llvm/Transforms/IPO/PassManagerBuilder.h"
static llvm::FunctionPassManager *TheFPM;
#endif

#define SUPPORT_TIMING 0
#if TIME
#define SUPPORT_TIMING 1
#include "copy_benchmark/hrtimer/hrtimer.h"
HRT_TIMESTAMP_T start, stop;
static uint64_t tmp;
static double commit_time = 0.0;
static double prof_time= 0.0;
#endif
#if SUPPORT_TIMING
unsigned long long g_timerfreq;
#endif

using namespace std;
using namespace llvm;

namespace farc {

static Module *TheModule;
static IRBuilder<> Builder(getGlobalContext());
static std::map<std::string, Value*> NamedValues;
static ExecutionEngine *TheExecutionEngine;
#define indent_str "  "

std::vector<std::string> Args;
FunctionType *FT;

#define LLVM_VOID     Type::getVoidTy(getGlobalContext())

#define LLVM_INT      Type::getInt32Ty(getGlobalContext())
#define LLVM_INT8     Type::getInt8Ty(getGlobalContext())
#define LLVM_INT32    Type::getInt32Ty(getGlobalContext())
#define LLVM_INT64    Type::getInt64Ty(getGlobalContext())
#define LLVM_INT8PTR  Type::getInt8PtrTy(getGlobalContext())


#define LLVM_FLOAT    Type::getFloatTy(getGlobalContext())
#define LLVM_DOUBLE   Type::getDoubleTy(getGlobalContext())

/** for debugging **
static Function* func; // printf function for debugging
std::vector<llvm::Type*> printf_arg_types;
printf_arg_types.push_back(LLVM_INT8PTR);
FunctionType* printf_type = FunctionType::get(LLVM_INT32, printf_arg_types, true);
Function *func = Function::Create(printf_type, Function::ExternalLinkage, Twine("printf"), TheModule);
Value *fmt_ptr = Builder.CreateGlobalStringPtr("stride to add: %i\n\0");
Value *fmt_ptr2 = Builder.CreateGlobalStringPtr("restore stride\n\0");
// now we can print as follows:
//llvm::CallInst *call = builder.CreateCall2(func, fmt_ptr, ValueToPrint);
*/

/* Codegen Functions */
static Value* multNode(int op1, Value* op2PtrNode) {
    Value* op1Node = ConstantInt::get(getGlobalContext(), APInt(64, op1, false));
    Value* op2Node = Builder.CreateIntCast(op2PtrNode, LLVM_INT64, false); 
    return Builder.CreateMul(op1Node, op2Node);
}

static inline Value* constNode(int val) {
    return ConstantInt::get(getGlobalContext(), APInt(32, val, false));
}

static inline Value* constNode(long val) {
    return ConstantInt::get(getGlobalContext(), APInt(64, val, false));
}

static inline void vectorCodegenUnrolled(Value* inbuf, Value* incount, Value* outbuf, Datatype* basetype, int count, int blocklen, int elemstride_in, int elemstride_out, bool pack) {
    Function *TheFunction = Builder.GetInsertBlock()->getParent();

    // Entry block
    Value* out = Builder.CreatePtrToInt(outbuf, LLVM_INT64);
    out->setName("out");
    Value* in = Builder.CreatePtrToInt(inbuf, LLVM_INT64);
    in->setName("in");


    // Outer loop
    BasicBlock *Preheader_outer_BB = Builder.GetInsertBlock();
    BasicBlock *Loop_outer_BB = BasicBlock::Create(getGlobalContext(), "outerloop", TheFunction);
    Builder.CreateBr(Loop_outer_BB);
    Builder.SetInsertPoint(Loop_outer_BB);

    // Induction var phi nodes
    PHINode *out1 = Builder.CreatePHI(LLVM_INT64, 2, "out1");
    out1->addIncoming(out, Preheader_outer_BB);
    PHINode *in1= Builder.CreatePHI(LLVM_INT64, 2, "in1");
    in1->addIncoming(in, Preheader_outer_BB);
    PHINode *i = Builder.CreatePHI(LLVM_INT32, 2, "i");
    i->addIncoming(constNode(0), Preheader_outer_BB);

    // Compute the size of the data written to the out buffer in the inner loop
    Value* nextin1 = NULL;
    Value* nextout1 = NULL;
    if (pack) {
        nextout1 = Builder.CreateAdd(out1, constNode(count * (long)elemstride_out));
	nextout1->setName("nextout1");
    } 
    else {
        nextin1 = Builder.CreateAdd(in1, constNode(count * (long)elemstride_in));
        nextin1->setName("nextin1");
    }

    // Inner loop
    BasicBlock *Preheader_inner_BB = Builder.GetInsertBlock();
    BasicBlock *Loop_inner_BB = BasicBlock::Create(getGlobalContext(), "innerloop", TheFunction);
    Builder.CreateBr(Loop_inner_BB);
    Builder.SetInsertPoint(Loop_inner_BB);
    
    // Induction var phi nodes
    PHINode *out2_1 = Builder.CreatePHI(LLVM_INT64, 2, "out2_1");
    out2_1->addIncoming(out1, Preheader_inner_BB);
    PHINode *in2_1 = Builder.CreatePHI(LLVM_INT64, 2, "in2_1");
    in2_1->addIncoming(in1, Preheader_inner_BB);
    
    // Cast out2 and in2 to pointers
    Value* out2_addr_1 = Builder.CreateIntToPtr(out2_1, LLVM_INT8PTR, "out2_addr_1");
    Value* in2_addr_1 = Builder.CreateIntToPtr(in2_1, LLVM_INT8PTR, "in2_addr_1");
    
    Value* out2_2 = Builder.CreateAdd(out2_1, constNode((long)elemstride_out), "out2_2");
    Value* in2_2 = Builder.CreateAdd(in2_1, constNode((long)elemstride_in), "in2_2");
    Value* out2_addr_2 = Builder.CreateIntToPtr(out2_2, LLVM_INT8PTR, "out2_addr_2");
    Value* in2_addr_2 = Builder.CreateIntToPtr(in2_2, LLVM_INT8PTR, "in2_addr_2");
    
    Value* out2_3 = Builder.CreateAdd(out2_2, constNode((long)elemstride_out), "out2_3");
    Value* in2_3 = Builder.CreateAdd(in2_2, constNode((long)elemstride_in), "in2_3");
    Value* out2_addr_3 = Builder.CreateIntToPtr(out2_3, LLVM_INT8PTR, "out2_addr_3");
    Value* in2_addr_3 = Builder.CreateIntToPtr(in2_3, LLVM_INT8PTR, "in2_addr_3");
    
    Value* out2_4 = Builder.CreateAdd(out2_3, constNode((long)elemstride_out), "out2_4");
    Value* in2_4 = Builder.CreateAdd(in2_3, constNode((long)elemstride_in), "in2_4");
    Value* out2_addr_4 = Builder.CreateIntToPtr(out2_4, LLVM_INT8PTR, "out2_addr_4");
    Value* in2_addr_4 = Builder.CreateIntToPtr(in2_4, LLVM_INT8PTR, "in2_addr_4");
    
    Value* out2_5 = Builder.CreateAdd(out2_4, constNode((long)elemstride_out), "out2_5");
    Value* in2_5 = Builder.CreateAdd(in2_4, constNode((long)elemstride_in), "in2_5");
    
    // Basetype Code Generation
    if (pack) {
        basetype->Codegen_Pack(in2_addr_1, ConstantInt::get(getGlobalContext(), APInt(32, blocklen, false)), out2_addr_1);
	basetype->Codegen_Pack(in2_addr_2, ConstantInt::get(getGlobalContext(), APInt(32, blocklen, false)), out2_addr_2);
	basetype->Codegen_Pack(in2_addr_3, ConstantInt::get(getGlobalContext(), APInt(32, blocklen, false)), out2_addr_3);
	basetype->Codegen_Pack(in2_addr_4, ConstantInt::get(getGlobalContext(), APInt(32, blocklen, false)), out2_addr_4);
    }
    else{
	basetype->Codegen_Unpack(in2_addr_1, ConstantInt::get(getGlobalContext(), APInt(32, blocklen, false)), out2_addr_1);
	basetype->Codegen_Unpack(in2_addr_2, ConstantInt::get(getGlobalContext(), APInt(32, blocklen, false)), out2_addr_2);
	basetype->Codegen_Unpack(in2_addr_3, ConstantInt::get(getGlobalContext(), APInt(32, blocklen, false)), out2_addr_3);
	basetype->Codegen_Unpack(in2_addr_4, ConstantInt::get(getGlobalContext(), APInt(32, blocklen, false)), out2_addr_4);
    }

    // check if we are finished with the loop over count
    Value* EndCond_inner = (pack) ? Builder.CreateICmpEQ(out2_5, nextout1, "innercond")
            : Builder.CreateICmpEQ(in2_5, nextin1, "innercond");

    // Create and branch to the inner loop postamble
    BasicBlock *LoopEnd_inner_BB = Builder.GetInsertBlock();
    BasicBlock *After_inner_BB = BasicBlock::Create(getGlobalContext(), "afterinner", TheFunction);
    Builder.CreateCondBr(EndCond_inner, After_inner_BB, Loop_inner_BB);
    Builder.SetInsertPoint(After_inner_BB);

    // Add backedges for the inner loop induction variables
    out2_1->addIncoming(out2_5, LoopEnd_inner_BB);
    in2_1->addIncoming(in2_5, LoopEnd_inner_BB);


    // Move the the extend-stride ptr back Extent(Basetype) * Stride - Size(Basetype) * Blocklen  
    if (pack) {
        nextin1 = Builder.CreateAdd(in1, constNode((long)(elemstride_in * (count-1) + elemstride_out)));
        nextin1->setName("nextin1");
    }
    else {
        nextout1 = Builder.CreateAdd(out1, constNode((long)(elemstride_out * (count-1) + elemstride_in)));
        nextout1->setName("nextout1");
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
    out1->addIncoming(nextout1, LoopEnd_outer_BB);
    in1->addIncoming(nextin1, LoopEnd_outer_BB);
    i->addIncoming(nexti, LoopEnd_outer_BB);    
}

static inline void vectorCodegenNormal(Value* inbuf, Value* incount, Value* outbuf, Datatype* basetype, int count, int blocklen, int elemstride_in, int elemstride_out, bool pack) {
    Function *TheFunction = Builder.GetInsertBlock()->getParent();

    // Entry block
    Value* out = Builder.CreatePtrToInt(outbuf, LLVM_INT64);
    out->setName("out");
    Value* in = Builder.CreatePtrToInt(inbuf, LLVM_INT64);
    in->setName("in");


    // Outer loop
    BasicBlock *Preheader_outer_BB = Builder.GetInsertBlock();
    BasicBlock *Loop_outer_BB = BasicBlock::Create(getGlobalContext(), "outerloop", TheFunction);
    Builder.CreateBr(Loop_outer_BB);
    Builder.SetInsertPoint(Loop_outer_BB);

    // Induction var phi nodes
    PHINode *out1 = Builder.CreatePHI(LLVM_INT64, 2, "out1");
    out1->addIncoming(out, Preheader_outer_BB);
    PHINode *in1= Builder.CreatePHI(LLVM_INT64, 2, "in1");
    in1->addIncoming(in, Preheader_outer_BB);
    PHINode *i = Builder.CreatePHI(LLVM_INT32, 2, "i");
    i->addIncoming(constNode(0), Preheader_outer_BB);
    
    // Compute the size of the data written to the out buffer in the inner loop
    Value* nextin1 = NULL;
    Value* nextout1 = NULL;
    if (pack) {
        nextout1 = Builder.CreateAdd(out1, constNode(count * (long)elemstride_out));
		nextout1->setName("nextout1");
    } 
    else {
        nextin1 = Builder.CreateAdd(in1, constNode(count * (long)elemstride_in));
		nextin1->setName("nextin1");
    }

    // Inner loop
    BasicBlock *Preheader_inner_BB = Builder.GetInsertBlock();
    BasicBlock *Loop_inner_BB = BasicBlock::Create(getGlobalContext(), "innerloop", TheFunction);
    Builder.CreateBr(Loop_inner_BB);
    Builder.SetInsertPoint(Loop_inner_BB);
    
    // Induction var phi nodes
    PHINode *out2 = Builder.CreatePHI(LLVM_INT64, 2, "out2");
    out2->addIncoming(out1, Preheader_inner_BB);
    PHINode *in2= Builder.CreatePHI(LLVM_INT64, 2, "in2");
    in2->addIncoming(in1, Preheader_inner_BB);
    
    // Cast out2 and in2 to pointers
    Value* out2_addr = Builder.CreateIntToPtr(out2, LLVM_INT8PTR);
    out2_addr->setName("out2_addr");
    Value* in2_addr = Builder.CreateIntToPtr(in2, LLVM_INT8PTR);
    in2_addr->setName("in2_addr");
    

    // Basetype Code Generation
    if (pack) basetype->Codegen_Pack(in2_addr, ConstantInt::get(getGlobalContext(), APInt(32, blocklen, false)), out2_addr);
    else      basetype->Codegen_Unpack(in2_addr, ConstantInt::get(getGlobalContext(), APInt(32, blocklen, false)), out2_addr);


    // Increment out2 and in2
    Value* nextout2 = Builder.CreateAdd(out2, constNode((long)elemstride_out));
    nextout2->setName("nextout2");
    Value* nextin2 = Builder.CreateAdd(in2, constNode((long)elemstride_in));
    nextin2->setName("nextin2");
    
    // check if we are finished with the loop over count
    Value* EndCond_inner = (pack) ? Builder.CreateICmpEQ(nextout2, nextout1, "innercond")
            : Builder.CreateICmpEQ(nextin2, nextin1, "innercond");
    
    // Create and branch to the inner loop postamble
    BasicBlock *LoopEnd_inner_BB = Builder.GetInsertBlock();
    BasicBlock *After_inner_BB = BasicBlock::Create(getGlobalContext(), "afterinner", TheFunction);
    Builder.CreateCondBr(EndCond_inner, After_inner_BB, Loop_inner_BB);
    Builder.SetInsertPoint(After_inner_BB);

    // Add backedges for the inner loop induction variables
    out2->addIncoming(nextout2, LoopEnd_inner_BB);
    in2->addIncoming(nextin2, LoopEnd_inner_BB);


    // Move the the extend-stride ptr back Extent(Basetype) * Stride - Size(Basetype) * Blocklen  
    if (pack) {
        nextin1 = Builder.CreateAdd(in1, constNode((long)(elemstride_in * (count-1) + elemstride_out)));
	nextin1->setName("nextin1");
    }
    else {
	nextout1 = Builder.CreateAdd(out1, constNode((long)(elemstride_out * (count-1) + elemstride_in)));
	nextout1->setName("nextout1");
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
    out1->addIncoming(nextout1, LoopEnd_outer_BB);
    in1->addIncoming(nextin1, LoopEnd_outer_BB);
    i->addIncoming(nexti, LoopEnd_outer_BB);
}

static void vectorCodegen(Value* inbuf, Value* incount, Value* outbuf, Datatype* basetype, int count, int blocklen, int elemstride_in, int elemstride_out, bool pack) {

    #if VECTOR_UNROLL
	if (count % 4 != 0) { 
    #endif
        vectorCodegenNormal(inbuf, incount, outbuf, basetype, count, blocklen, elemstride_in, elemstride_out, pack);
    #if VECTOR_UNROLL
    }
    else {
        vectorCodegenUnrolled(inbuf, incount, outbuf, basetype, count, blocklen, elemstride_in, elemstride_out, pack);
    }
    #endif
}


/* Datatype */
void Datatype::print() {
	printf("%s\n", this->toString().c_str());
}


/* PrimitiveDatatype */
PrimitiveDatatype::PrimitiveDatatype(PrimitiveDatatype::PrimitiveType type) : Datatype() {

    this->Type = type;

    if (Type == BYTE)   this->Extent = 1;
    if (Type == CHAR)   this->Extent = 1;
    if (Type == DOUBLE) this->Extent = sizeof(double);
    if (Type == FLOAT)  this->Extent = sizeof(float);
    if (Type == INT)    this->Extent = sizeof(int);
    //TODO add more. Remember to also add them to the print function.

    this->Size = this->Extent;

}

PrimitiveDatatype* PrimitiveDatatype::Clone() {

    PrimitiveDatatype* t_new = new PrimitiveDatatype(this->Type);
    return t_new;
}

static inline void vmove(Value *dst, Value *src, int count, llvm::Type *elemtype) {
	llvm::Type *elemvectype_ptr = PointerType::getUnqual(VectorType::get(elemtype, count));
	Value *in_vec = Builder.CreateBitCast(src, elemvectype_ptr, "in2_addr_vec");
	Value *out_vec = Builder.CreateBitCast(dst, elemvectype_ptr, "out2_addr_vec");
	Value *elems = Builder.CreateAlignedLoad(in_vec, 1, "elems");
	Builder.CreateAlignedStore(elems, out_vec, 1);
}

static inline Value *incrementPtr(Value *ptr, int byteInc) {
	Value *addr = Builder.CreatePtrToInt(ptr, LLVM_INT64);
	Value *newaddr = Builder.CreateAdd(addr, Builder.getInt64(byteInc));
	return Builder.CreateIntToPtr(newaddr, LLVM_INT8PTR);
}

static inline llvm::Type *toLLVMType(PrimitiveDatatype::PrimitiveType type) {
	llvm::Type *elemtype = NULL;
	switch (type) {
	case PrimitiveDatatype::DOUBLE:
		elemtype = LLVM_DOUBLE;
		break;
	case PrimitiveDatatype::INT:
		elemtype = LLVM_INT;
		break;
	case PrimitiveDatatype::FLOAT:
		elemtype = LLVM_FLOAT;
		break;
	case PrimitiveDatatype::BYTE:
		elemtype = LLVM_INT8;
		break;
	case PrimitiveDatatype::CHAR:
		elemtype = LLVM_INT8;
		break;
	default:
		fprintf(stderr, "Type not supported");
		assert(false);
	}
	assert(elemtype != NULL);
	return elemtype;
}

void PrimitiveDatatype::Codegen_Pack(Value* inbuf, Value* incount, Value* outbuf) {

    Function* TheFunction = Builder.GetInsertBlock()->getParent();
    llvm::ConstantInt* incount_ci = dyn_cast<llvm::ConstantInt>(incount);
    // Case 3
    if (incount_ci != NULL) {

#if !PACKVAR
#undef PACKVAR
#define PACKVAR 8
#endif

// NOTE that most of these packvars don't work with unaligned memory
// they are only here for benchmarking/comparing different variants
#if PACKVAR == 1        
        // Bitcast instructions that make it easier to interface with the outside code.
        // Note that these don't result in any assembly instructions
        llvm::Type* vectypeptr = PointerType::getUnqual(VectorType::get(LLVM_INT8, this->getSize() * incount_ci->getSExtValue()));
        Value* out_vec = Builder.CreateBitCast(outbuf, vectypeptr, "out2_addr_vec");
        Value* in_vec = Builder.CreateBitCast(inbuf, vectypeptr, "in2_addr_vec");

        // Load and store
        Value* bytes = Builder.CreateAlignedLoad(in_vec, 1, "bytes");
        Builder.CreateAlignedStore(bytes, out_vec, 1);
#elif PACKVAR == 2       
        // Bitcast instructions that make it easier to interface with the outside code.
        // Note that these don't result in any assembly instructions
        llvm::Type* vectypeptr = PointerType::getUnqual(VectorType::get(LLVM_INT8, this->getSize() * incount_ci->getSExtValue()));
        Value* out_vec = Builder.CreateBitCast(outbuf, vectypeptr, "out2_addr_vec");
        Value* in_vec = Builder.CreateBitCast(inbuf, vectypeptr, "in2_addr_vec");

        // Load and store
        Value* bytes = Builder.CreateAlignedLoad(in_vec, 16, "bytes");
        Builder.CreateAlignedStore(bytes, out_vec, 16);
#elif PACKVAR == 3
        // Bitcast instructions that make it easier to interface with the outside code.
        // Note that these don't result in any assembly instructions
        llvm::Type* vectypeptr = PointerType::getUnqual(VectorType::get(LLVM_INT8, this->getSize() * incount_ci->getSExtValue()));
        Value* out_vec = Builder.CreateBitCast(outbuf, vectypeptr, "out2_addr_vec");
        Value* in_vec = Builder.CreateBitCast(inbuf, vectypeptr, "in2_addr_vec");

        // Load and store
		MDNode *Node = MDNode::get(getGlobalContext(), Builder.getInt32(1));
        Value* bytes = Builder.CreateAlignedLoad(in_vec, 16, "bytes");
        StoreInst* store = Builder.CreateAlignedStore(bytes, out_vec, 16);
		store->setMetadata(TheModule->getMDKindID("nontemporal"), Node); // why the hell do i have to supply 1 here llvm?
#elif PACKVAR == 4
		int size_to_pack = this->getSize() * incount_ci->getSExtValue();
        llvm::Type* vectypeptr = PointerType::getUnqual(VectorType::get(Type::getDoubleTy(getGlobalContext()), 2));

		while (size_to_pack > 0) {
		
			MDNode *Node = MDNode::get(getGlobalContext(), Builder.getInt32(1));
			Value* out_vec = Builder.CreateBitCast(outbuf, vectypeptr, "out2_addr_vec");
			Value* in_vec = Builder.CreateBitCast(inbuf, vectypeptr, "in2_addr_vec");
			Value* bytes = Builder.CreateAlignedLoad(in_vec, 16, "bytes");
			StoreInst* store = Builder.CreateAlignedStore(bytes, out_vec, 16);
			store->setMetadata(TheModule->getMDKindID("nontemporal"), Node); // why the hell do i have to supply 1 here llvm?

    		Value* out_addr_cvi = Builder.CreatePtrToInt(outbuf, LLVM_INT64);
    		Value* out_addr = Builder.CreateAdd(out_addr_cvi,  Builder.getInt64(16));
    		outbuf = Builder.CreateIntToPtr(out_addr, LLVM_INT8PTR);

    		Value* in_addr_cvi = Builder.CreatePtrToInt(inbuf, LLVM_INT64);
    		Value* in_addr = Builder.CreateAdd(in_addr_cvi, Builder.getInt64(16));
    		inbuf = Builder.CreateIntToPtr(in_addr, LLVM_INT8PTR);
			
			size_to_pack -= 16;
		}

		// revert changes to inbuf, outbuf
    	Value* out_addr_cvi = Builder.CreatePtrToInt(outbuf, LLVM_INT64);
    	Value* out_addr = Builder.CreateAdd(out_addr_cvi, Builder.getInt64(-size_to_pack));
    	outbuf = Builder.CreateIntToPtr(out_addr, LLVM_INT8PTR);

    	Value* in_addr_cvi = Builder.CreatePtrToInt(inbuf, LLVM_INT64);
    	Value* in_addr = Builder.CreateAdd(in_addr_cvi, Builder.getInt64(-size_to_pack));
    	inbuf = Builder.CreateIntToPtr(in_addr, LLVM_INT8PTR);
#elif PACKVAR == 5
		int size_to_pack = this->getSize() * incount_ci->getSExtValue();
        llvm::Type* vectypeptr = PointerType::getUnqual(VectorType::get(Type::getDoubleTy(getGlobalContext()), 2));

		while (size_to_pack > 0) {
		
			MDNode *Node = MDNode::get(getGlobalContext(), Builder.getInt32(1));
			Value* out_vec = Builder.CreateBitCast(outbuf, vectypeptr, "out2_addr_vec");
			Value* in_vec = Builder.CreateBitCast(inbuf, vectypeptr, "in2_addr_vec");
			Value* bytes = Builder.CreateAlignedLoad(in_vec, 16, "bytes");
			StoreInst* store = Builder.CreateAlignedStore(bytes, out_vec, 16);

    		Value* out_addr_cvi = Builder.CreatePtrToInt(outbuf, LLVM_INT64);
    		Value* out_addr = Builder.CreateAdd(out_addr_cvi,  Builder.getInt64(16));
    		outbuf = Builder.CreateIntToPtr(out_addr, LLVM_INT8PTR);

    		Value* in_addr_cvi = Builder.CreatePtrToInt(inbuf, LLVM_INT64);
    		Value* in_addr = Builder.CreateAdd(in_addr_cvi, Builder.getInt64(16));
    		inbuf = Builder.CreateIntToPtr(in_addr, LLVM_INT8PTR);
			
			size_to_pack -= 16;
		}

		// revert changes to inbuf, outbuf
    	Value* out_addr_cvi = Builder.CreatePtrToInt(outbuf, LLVM_INT64);
    	Value* out_addr = Builder.CreateAdd(out_addr_cvi, Builder.getInt64(-size_to_pack));
    	outbuf = Builder.CreateIntToPtr(out_addr, LLVM_INT8PTR);

    	Value* in_addr_cvi = Builder.CreatePtrToInt(inbuf, LLVM_INT64);
    	Value* in_addr = Builder.CreateAdd(in_addr_cvi, Builder.getInt64(-size_to_pack));
    	inbuf = Builder.CreateIntToPtr(in_addr, LLVM_INT8PTR);
#elif PACKVAR == 6
		Value* contig_extend = multNode(this->getSize(), incount);
        Value* memcopy = Builder.CreateMemCpy(outbuf, inbuf, contig_extend, 1);
#elif PACKVAR == 7
        int size_to_pack = this->getSize() * incount_ci->getSExtValue();

        // if we know (at compile time) there are only few bytes (not enough for vector insts), just codegen/unroll unaligned instructions to copy them
        if ((size_to_pack < 16) || ((size_to_pack < 32) && (size_to_pack & ~15))) {
            while (size_to_pack > 0) {
                int pack_now = 0;
                if (size_to_pack > 7) pack_now = 8;
                else if (size_to_pack > 3) pack_now = 4;
                else pack_now = 1;

                // copy packed bytes
                llvm::Type* vectypeptr = PointerType::getUnqual(VectorType::get(Type::getInt8Ty(getGlobalContext()), pack_now));
                Value* out_vec = Builder.CreateBitCast(outbuf, vectypeptr, "smallcopy_out_vec");
                Value* in_vec = Builder.CreateBitCast(inbuf, vectypeptr, "smallcopy_in_vec");
                Value* bytes = Builder.CreateLoad(in_vec, "bytes");
                Builder.CreateStore(bytes, out_vec);

                //increment inbuf and outbuf by "packed"
                Value* in_addr_cvi = Builder.CreatePtrToInt(inbuf, LLVM_INT64);
                Value* in_addr = Builder.CreateAdd(in_addr_cvi, Builder.getInt64(pack_now));
                inbuf = Builder.CreateIntToPtr(in_addr, LLVM_INT8PTR);
                Value* out_addr_cvi = Builder.CreatePtrToInt(outbuf, LLVM_INT64);
                Value* out_addr = Builder.CreateAdd(out_addr_cvi,  Builder.getInt64(pack_now));
                outbuf = Builder.CreateIntToPtr(out_addr, LLVM_INT8PTR);

                size_to_pack -= pack_now;
            }
        }
        else {

            // copy single bytes until outbuf is 16 byte aligned
            // this has to be a test-first loop, since outbuf could already be aligned

            Value* size_to_pack = constNode((long) (this->getSize() * incount_ci->getSExtValue()));

            Value* in = Builder.CreatePtrToInt(inbuf, LLVM_INT64);
            Value* out = Builder.CreatePtrToInt(outbuf, LLVM_INT64);

            // loop
            BasicBlock *Preheader_prefix_BB = Builder.GetInsertBlock();
            BasicBlock *Condition_prefix_BB = BasicBlock::Create(getGlobalContext(), "prefixcondition", TheFunction);
            BasicBlock *Loop_prefix_BB = BasicBlock::Create(getGlobalContext(), "prefixloop", TheFunction);
            BasicBlock *After_prefix_BB = BasicBlock::Create(getGlobalContext(), "afterprefix", TheFunction);
            Builder.CreateBr(Condition_prefix_BB);

            Builder.SetInsertPoint(Condition_prefix_BB);

            // Induction var phi nodes
            PHINode *out2 = Builder.CreatePHI(LLVM_INT64, 2, "out2");
            out2->addIncoming(out, Preheader_prefix_BB);
            PHINode *in2= Builder.CreatePHI(LLVM_INT64, 2, "in2");
            in2->addIncoming(in, Preheader_prefix_BB);
            PHINode *size_to_pack2 = Builder.CreatePHI(LLVM_INT64, 2, "size_to_pack2");
            size_to_pack2->addIncoming(size_to_pack, Preheader_prefix_BB);

            Value* out_masked = Builder.CreateAnd(out2, constNode(0xFL));
            Value* StartCond_prefix = Builder.CreateICmpEQ(out_masked, constNode(0L), "prefixstartcond");
            Builder.CreateCondBr(StartCond_prefix, After_prefix_BB, Loop_prefix_BB);
            Builder.SetInsertPoint(Loop_prefix_BB);

            // Cast out2 and in2 to pointers
            Value* out2_addr = Builder.CreateIntToPtr(out2, LLVM_INT8PTR, "out2_addr");
            Value* in2_addr = Builder.CreateIntToPtr(in2, LLVM_INT8PTR, "in2_addr");

            //load-store
            Value* byte = Builder.CreateLoad(in2_addr, "byte");
            Builder.CreateStore(byte, out2_addr);

            // Increment out2 and in2, decrement next_size_to_pack
            Value* nextout2 = Builder.CreateAdd(out2, constNode(1L), "nextout2");
            Value* nextin2 = Builder.CreateAdd(in2, constNode(1L), "nextin2");
            Value* next_size_to_pack = Builder.CreateSub(size_to_pack2, constNode(1L), "next_size_to_pack");

            // Create and branch to the prefix loop postamble
            Builder.CreateBr(Condition_prefix_BB);

            // Add backedges for the prefix loop induction variables
            out2->addIncoming(nextout2, Loop_prefix_BB);
            in2->addIncoming(nextin2, Loop_prefix_BB);
            size_to_pack2->addIncoming(next_size_to_pack, Loop_prefix_BB);
            Builder.SetInsertPoint(After_prefix_BB);

           
            // outbuf is now 16byte aligned, copy as much as possible with aligned stores

            // loop
            BasicBlock *Preheader_aligned_BB = Builder.GetInsertBlock();
            BasicBlock *Condition_aligned_BB = BasicBlock::Create(getGlobalContext(), "alignedcondition", TheFunction);
            BasicBlock *Loop_aligned_BB = BasicBlock::Create(getGlobalContext(), "alignedloop", TheFunction);
            BasicBlock *After_aligned_BB = BasicBlock::Create(getGlobalContext(), "afteraligned", TheFunction);
            Builder.CreateBr(Condition_aligned_BB);

            Builder.SetInsertPoint(Condition_aligned_BB);

            // Induction var phi nodes
            PHINode *out_aligned_2 = Builder.CreatePHI(LLVM_INT64, 2, "out_aligned_2");
            out_aligned_2->addIncoming(out2, Preheader_aligned_BB);
            PHINode *in_aligned_2= Builder.CreatePHI(LLVM_INT64, 2, "in_aligned_2");
            in_aligned_2->addIncoming(in2, Preheader_aligned_BB);
            PHINode *size_to_pack_aligned_2 = Builder.CreatePHI(LLVM_INT64, 2, "size_to_pack_aligned_2");
            size_to_pack_aligned_2->addIncoming(size_to_pack2, Preheader_aligned_BB);

            Value* StartCond_aligned = Builder.CreateICmpULT(size_to_pack_aligned_2, constNode(16L), "alignedstartcond");
            Builder.CreateCondBr(StartCond_aligned, After_aligned_BB, Loop_aligned_BB);
            Builder.SetInsertPoint(Loop_aligned_BB);

            // Cast out_aligned_2 and in_aligned_2 to pointers
            Value* out_aligned_2_addr = Builder.CreateIntToPtr(out_aligned_2, LLVM_INT8PTR, "out_aligned_2_addr");
            Value* in_aligned_2_addr = Builder.CreateIntToPtr(in_aligned_2, LLVM_INT8PTR, "in_aligned_2_addr");

            //load-store
            llvm::Type* vectypeptr = PointerType::getUnqual(VectorType::get(Type::getDoubleTy(getGlobalContext()), 2));
            Value* out_aligned_vec = Builder.CreateBitCast(out_aligned_2_addr, vectypeptr, "out2_addr_vec");
			Value* in_aligned_vec = Builder.CreateBitCast(in_aligned_2_addr, vectypeptr, "in2_addr_vec");
			Value* bytes_aligned = Builder.CreateLoad(in_aligned_vec, "bytes_aligned");
			Builder.CreateAlignedStore(bytes_aligned, out_aligned_vec, 16);

            // Increment out_aligned_2 and in_aligned_2, decrement size_to_pack_aligned_2
            Value* nextout_aligned_2 = Builder.CreateAdd(out_aligned_2, constNode(16L), "nextout_aligned_2");
            Value* nextin_aligned_2 = Builder.CreateAdd(in_aligned_2, constNode(16L), "nextin_aligned_2");
            Value* next_size_to_pack_aligned = Builder.CreateSub(size_to_pack_aligned_2, constNode(16L), "next_size_to_pack_aligned");

            // Create and branch to the aligned loop postamble
            Builder.CreateBr(Condition_aligned_BB);

            // Add backedges for the aligned loop induction variables
            out_aligned_2->addIncoming(nextout_aligned_2, Loop_aligned_BB);
            in_aligned_2->addIncoming(nextin_aligned_2, Loop_aligned_BB);
            size_to_pack_aligned_2->addIncoming(next_size_to_pack_aligned, Loop_aligned_BB);
            Builder.SetInsertPoint(After_aligned_BB);


            // copy the remaining bytes in an unaligned manner

            BasicBlock *Preheader_tail_BB = Builder.GetInsertBlock();
            BasicBlock *Condition_tail_BB = BasicBlock::Create(getGlobalContext(), "tailcondition", TheFunction);
            BasicBlock *Loop_tail_BB = BasicBlock::Create(getGlobalContext(), "tailloop", TheFunction);
            BasicBlock *After_tail_BB = BasicBlock::Create(getGlobalContext(), "aftertail", TheFunction);
            Builder.CreateBr(Condition_tail_BB);

            Builder.SetInsertPoint(Condition_tail_BB);

            // Induction var phi nodes
            PHINode *out_tail_2 = Builder.CreatePHI(LLVM_INT64, 2, "out_tail_2");
            out_tail_2->addIncoming(out_aligned_2, Preheader_tail_BB);
            PHINode *in_tail_2= Builder.CreatePHI(LLVM_INT64, 2, "in_tail_2");
            in_tail_2->addIncoming(in_aligned_2, Preheader_tail_BB);
            PHINode *size_to_pack_tail_2 = Builder.CreatePHI(LLVM_INT64, 2, "size_to_pack_tail_2");
            size_to_pack_tail_2->addIncoming(size_to_pack_aligned_2, Preheader_tail_BB);

            Value* StartCond_tail = Builder.CreateICmpEQ(size_to_pack_tail_2, constNode(0L), "tailstartcond");
            Builder.CreateCondBr(StartCond_tail, After_tail_BB, Loop_tail_BB);
            Builder.SetInsertPoint(Loop_tail_BB);

            // Cast out_tail_2 and in_tail_2 to pointers
            Value* out_tail_2_addr = Builder.CreateIntToPtr(out_tail_2, LLVM_INT8PTR, "out_tail_2_addr");
            Value* in_tail_2_addr = Builder.CreateIntToPtr(in_tail_2, LLVM_INT8PTR, "in_tail_2_addr");

            //load-store
            Value* byte_tail = Builder.CreateLoad(in_tail_2_addr, "byte");
            Builder.CreateStore(byte_tail, out_tail_2_addr);

            // Increment out_tail_2 and in_tail_2, decrement size_to_pack_tail_2
            Value* nextout_tail_2 = Builder.CreateAdd(out_tail_2, constNode(1L), "nextout_tail_2");
            Value* nextin_tail_2 = Builder.CreateAdd(in_tail_2, constNode(1L), "nextin_tail_2");
            Value* next_size_to_pack_tail = Builder.CreateSub(size_to_pack_tail_2, constNode(1L), "next_size_to_pack_tail");

            // Create and branch to the tail loop postamble
            Builder.CreateBr(Condition_tail_BB);

            // Add backedges for the tail loop induction variables
            out_tail_2->addIncoming(nextout_tail_2, Loop_tail_BB);
            in_tail_2->addIncoming(nextin_tail_2, Loop_tail_BB);
            size_to_pack_tail_2->addIncoming(next_size_to_pack_tail, Loop_tail_BB);
            Builder.SetInsertPoint(After_tail_BB);

        }
#elif PACKVAR == 8
		// Kernel that performs unaligned memcopies using vector
		// instructions.  The kernel copies as many elements as
		// possible using vector instructions of a size specified by
		// SIMD_BYTE_SIZE in a loop that is unrolled COPY_LOOP_UNROLL
		// times.
		//
		// If the size of the elements is not divisible by
		// SIMD_BYTE_SIZE then the overflow elements are copied using
		// succesively smaller pow-of-two sized vector instructions in
		// a postamble until all are copied.
		//
		// If the number of bytes to copy is less than
		// COPY_LOOP_TRESHOLD then the code skips the loop creating
		// all-together and falls through to the postamble generator,
		// which produces fully unrolled code.

		// This should be a power of two (on kepler performance is bad
		// with 32-byte simd)
		#ifndef SIMD_BYTE_SIZE
		#define SIMD_BYTE_SIZE 16
		#endif

		// Number of bytes to allow before introducing a loop.
		// Tip: Make this 4*8 times larger than COPY_LOOP_UNROLL for
		// smooth double performance
		#ifndef COPY_LOOP_TRESHOLD
		#define COPY_LOOP_TRESHOLD 32*8
		#endif

		// Number of times to unroll the copy loop
		#ifndef COPY_LOOP_UNROLL
		#define COPY_LOOP_UNROLL 8
		#endif

		const int LOOP_ELEM_TRESHOLD = COPY_LOOP_TRESHOLD / this->Size;
		const int vector_size = SIMD_BYTE_SIZE / this->Size;

		assert((SIMD_BYTE_SIZE & (SIMD_BYTE_SIZE-1)) == 0); // Assert power-of-two
		assert((SIMD_BYTE_SIZE % this->Size) == 0);
		assert((vector_size & (vector_size-1)) == 0); // Assert power-of-two

		llvm::Type *elemtype = toLLVMType(this->Type);

		// Number of elements to copy
		int incount_val = incount_ci->getSExtValue();

		// Trunk to a multiple of the unroll factor
		const int vectors_to_copy = ((incount_val / vector_size) / COPY_LOOP_UNROLL) * COPY_LOOP_UNROLL;

		// Copy vectors_to_copy  vectors of size vecsize.
		// If the number of values is less than the LOOP_ELEM_TRESHOLD
		// treshold then we fall through to the postamble and unroll
		// everything
		if (vectors_to_copy > 0 && incount_val >= LOOP_ELEM_TRESHOLD) {
			Value *inbuf_int = Builder.CreatePtrToInt(inbuf, LLVM_INT64);
			Value *exitval   = Builder.CreateAdd(inbuf_int, constNode((long)vectors_to_copy * SIMD_BYTE_SIZE), "exitval");

			BasicBlock *header = Builder.GetInsertBlock();
			BasicBlock *copyloop   = BasicBlock::Create(getGlobalContext(), "copyloop", TheFunction);
			Builder.CreateBr(copyloop);
			Builder.SetInsertPoint(copyloop);

			PHINode *inphi = Builder.CreatePHI(LLVM_INT8PTR, 2, "in3");
			inphi->addIncoming(inbuf, header);
			PHINode *outphi  = Builder.CreatePHI(LLVM_INT8PTR, 2, "out3");
			outphi->addIncoming(outbuf, header);
			inbuf  = inphi;
			outbuf = outphi;

			Value *in_addr = NULL;
			for (int i=0; i<COPY_LOOP_UNROLL; i++) {
				vmove(outbuf, inbuf, vector_size, elemtype);

				Value *in_addr_cvi = Builder.CreatePtrToInt(inbuf, LLVM_INT64);
				in_addr = Builder.CreateAdd(in_addr_cvi, Builder.getInt64(this->Size * vector_size));
				inbuf = Builder.CreateIntToPtr(in_addr, LLVM_INT8PTR);

				outbuf = incrementPtr(outbuf, this->Size * vector_size);
			}
			incount_val -= vectors_to_copy * vector_size;

			inphi->addIncoming(inbuf, copyloop);
			outphi->addIncoming(outbuf, copyloop);

			// Create and jump to postamble
			BasicBlock *copypostamble = BasicBlock::Create(getGlobalContext(), "copypostamble", TheFunction);
			Value *exitcond = Builder.CreateICmpEQ(in_addr, exitval);
			Builder.CreateCondBr(exitcond, copypostamble, copyloop);
			Builder.SetInsertPoint(copypostamble);
		}

		// Copy postamble: copy the overflow elements that did not fit in full vector
		for (int vecsize=vector_size; vecsize > 0; vecsize /= 2) {
			const int veccount = incount_val / vecsize;
			for (int i=0; i<veccount; i++) {
				vmove(outbuf, inbuf, vecsize, elemtype);
				inbuf  = incrementPtr(inbuf, this->Size * vecsize);
				outbuf = incrementPtr(outbuf, this->Size * vecsize);
			}
			incount_val -= veccount * vecsize;
		}
#else
#error NO PACKVAR DEFINED
#endif

    }
    // Case 1 and 2
    else {
		Value* contig_extend = multNode(this->getSize(), incount);
        Builder.CreateMemCpy(outbuf, inbuf, contig_extend, 1);
    }

}

void PrimitiveDatatype::Codegen_Unpack(Value* inbuf, Value* incount, Value* outbuf) {
    // does exactly the same as pack for primitive types
    return this->Codegen_Pack(inbuf, incount, outbuf);
}

int PrimitiveDatatype::getExtent() {
    return this->Extent;
}

int PrimitiveDatatype::getSize() {
    return this->Size;
}

string PrimitiveDatatype::toString() {
	string res = "";
    switch (this->Type) {
    case BYTE:
        res.append("byte");
        break;
    case CHAR:
        res.append("char");
        break;
    case DOUBLE:
        res.append("double");
        break;
    case FLOAT:
        res.append("float");
        break;
    case INT:
        res.append("int");
        break;
    default:
        res.append("N/A");
        break;
    }
	return res;
}

/* Value* PrimitiveDatatype::Codegen_Pack_partial(Value* inbuf, Value* incount, Value* outbuf, Value* outbuf_from, Value* outbuf_to) {

    // TODO pseudocode for this function 

    //TODO we shouldn't "cut" primitive types, but this pseudo-code
    //can still be used for other ddts

    if (incount * this->size() < outbuf_from) {
        // case 1: we don't start packing at this node
        inbuf += incount * this->getExtent();
        outbuf += incount * this->getSize();
        return;
    }

    // jump over the blocks which are out of range
    x = outbuf_from / this->getSize();
    inbuf += x * this->getExtent();
    outbuf += x * this->getSize();
    if (outbuf != outbuf_from) {
        // pack a fraction of a primitive type 
        // (easy in this case, since we know it's contiguous)
        trunc_front = outbuf_from - outbuf;
        trunc_back = this->getSize() - trunc_front;
        memcpy(outbuf, inbuf, trunc_back);
        inbuf += trunc_back;
        outbuf += trunc_back;
        incount -= x + 1;
    }

    // generate fast code for as many as possible "normal" blocks
    //TODO calc fullblocks
    fullblocks = (outbuf_to - outbuf_from) / 
    Codegen_Pack(inbuf, fullblocks, outbuf);
    inbuf += fullblocks * this->getExtent();
    outbuf += fullblocks * this->getSize();
    incount -= x + fullblocks;

    // return how many bytes have been packed 

} */


/* Class ContiguousDatatype */
ContiguousDatatype::ContiguousDatatype(Datatype* type, int count) {

    this->Basetype = type->Clone();
    this->Count = count;

}

ContiguousDatatype::~ContiguousDatatype(void) {
    
    delete(this->Basetype);

}

ContiguousDatatype* ContiguousDatatype::Clone() {

    ContiguousDatatype* t_new = new ContiguousDatatype(this->Basetype, this->Count);

    return t_new;

}

void ContiguousDatatype::Codegen(Value* inbuf, Value* incount, Value* outbuf, int elemstride_in, int elemstride_out, bool pack) {
    Function* TheFunction = Builder.GetInsertBlock()->getParent();

    // Loop
    BasicBlock* PreheaderBB = Builder.GetInsertBlock();
    BasicBlock* LoopBB = BasicBlock::Create(getGlobalContext(), "loop", TheFunction);
    Builder.CreateBr(LoopBB);
    Builder.SetInsertPoint(LoopBB);

    // Induction var phi nodes
    PHINode *out = Builder.CreatePHI(LLVM_INT8PTR, 2, "out");
    out->addIncoming(outbuf, PreheaderBB);
    PHINode *in= Builder.CreatePHI(LLVM_INT8PTR, 2, "in");
    in->addIncoming(inbuf, PreheaderBB);
    PHINode *i = Builder.CreatePHI(LLVM_INT32, 2, "i");
    i->addIncoming(constNode(0), PreheaderBB);


    // Basetype Code Generation
    if (pack) Basetype->Codegen_Pack(in, ConstantInt::get(getGlobalContext(), APInt(32, this->Count, false)), out);
    else      Basetype->Codegen_Unpack(in, ConstantInt::get(getGlobalContext(), APInt(32, this->Count, false)), out);


    // Increment the out ptr by Size(Basetype) * Blocklen
    Value* out_bytes_to_stride = ConstantInt::get(getGlobalContext(), APInt(64, elemstride_out * this->Count, false));
    Value* out_addr_cvi = Builder.CreatePtrToInt(out, LLVM_INT64);
    Value* out_addr = Builder.CreateAdd(out_addr_cvi, out_bytes_to_stride);
    Value* nextout = Builder.CreateIntToPtr(out_addr, LLVM_INT8PTR);

    // Increment the in ptr by Extent(Basetype) * Stride
    Value* in_bytes_to_stride = ConstantInt::get(getGlobalContext(), APInt(64, elemstride_in * this->Count, false));
    Value* in_addr_cvi = Builder.CreatePtrToInt(in, LLVM_INT64);
    Value* in_addr = Builder.CreateAdd(in_addr_cvi, in_bytes_to_stride);
    Value* nextin = Builder.CreateIntToPtr(in_addr, LLVM_INT8PTR);

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

void ContiguousDatatype::Codegen_Pack(Value* inbuf, Value* incount, Value* outbuf) {
    Codegen(inbuf, incount, outbuf, this->Basetype->getExtent(), this->Basetype->getSize(), true);
}

void ContiguousDatatype::Codegen_Unpack(Value* inbuf, Value* incount, Value* outbuf) {
    Codegen(inbuf, incount, outbuf, this->Basetype->getSize(), this->Basetype->getExtent(), false);
}

int ContiguousDatatype::getExtent() {
    return this->Count * this->Basetype->getExtent();
}

int ContiguousDatatype::getSize() {
    return this->Count * this->Basetype->getSize();
}

string ContiguousDatatype::toString() {
	stringstream res;
	res << "ctg(" << this->Count << ")[" << Basetype->toString() << "]";
	return res.str();
}


/* Class VectorDatatype */
VectorDatatype::VectorDatatype(Datatype* type, int count, int blocklen, int stride) {

    this->Basetype = type->Clone();
    this->Count = count;
    this->Blocklen = blocklen;
    this->Stride = stride;

} 

VectorDatatype::~VectorDatatype(void) {

    delete(this->Basetype);

}

VectorDatatype* VectorDatatype::Clone() {

    VectorDatatype* t_new = new VectorDatatype(this->Basetype, this->Count, this->Blocklen, this->Stride);

    return t_new;

}

void VectorDatatype::Codegen_Pack(Value* inbuf, Value* incount, Value* outbuf) {
    vectorCodegen(inbuf, incount, outbuf, this->Basetype, this->Count,
            this->Blocklen, this->Basetype->getExtent() * this->Stride, this->Basetype->getSize() * this->Blocklen, true);
}

void VectorDatatype::Codegen_Unpack(Value* inbuf, Value* incount, Value* outbuf) {
    vectorCodegen(inbuf, incount, outbuf, this->Basetype, this->Count, this->Blocklen,
            this->Basetype->getSize() * this->Blocklen, this->Basetype->getExtent() * this->Stride, false);
}

int VectorDatatype::getExtent() {
    return (this->Count - 1) * this->Basetype->getExtent() * this->Stride + this->Blocklen * this->Basetype->getExtent();
}

int VectorDatatype::getSize() {
    return this->Count * this->Blocklen*this->Basetype->getSize();
}

string VectorDatatype::toString() {
	stringstream res;
	res << "vec(" << this->Count << " " << this->Blocklen << " " << this->Stride << ")[" << Basetype->toString() << "]";
	return res.str();
}


/* Class HVectorDatatype */
HVectorDatatype::HVectorDatatype(Datatype* type, int count, int blocklen, int stride) {

    this->Basetype = type->Clone();
    this->Count = count;
    this->Blocklen = blocklen;
    this->Stride = stride;

} 

HVectorDatatype* HVectorDatatype::Clone() {

    HVectorDatatype* t_new = new HVectorDatatype(this->Basetype, this->Count, this->Blocklen, this->Stride);

    return t_new;

}

HVectorDatatype::~HVectorDatatype(void) {

    delete(this->Basetype);

}

void HVectorDatatype::Codegen_Pack(Value* inbuf, Value* incount, Value* outbuf) {
    vectorCodegen(inbuf, incount, outbuf, this->Basetype, this->Count, this->Blocklen,
            this->Stride, this->Basetype->getSize() * this->Blocklen, true);
}

void HVectorDatatype::Codegen_Unpack(Value* inbuf, Value* incount, Value* outbuf) {
    vectorCodegen(inbuf, incount, outbuf, this->Basetype, this->Count, this->Blocklen,
            this->Basetype->getSize() * this->Blocklen, this->Stride, false);
}

int HVectorDatatype::getExtent() {
    return (this->Count-1)*this->Stride + this->Blocklen*this->Basetype->getExtent();
}

int HVectorDatatype::getSize() {
    return this->Count * this->Blocklen*this->Basetype->getSize();
}

string HVectorDatatype::toString() {
	stringstream res;
	res << "hvec(" << this->Count << " " << this->Blocklen << " " << this->Stride << ")[" << Basetype->toString() << "]";
	return res.str();
}


/* Class IndexedBlockDatatype */
IndexedBlockDatatype::IndexedBlockDatatype(int count, int blocklen, int* displ, Datatype* basetype) : Datatype() {

    this->Count = count;
    this->Basetype = basetype->Clone();
    this->Blocklen = blocklen;
    for (int i=0; i<count; i++) this->Displ.push_back(displ[i]);

}

IndexedBlockDatatype::~IndexedBlockDatatype(void) {

    delete(this->Basetype);

}

IndexedBlockDatatype* IndexedBlockDatatype::Clone() {

    IndexedBlockDatatype* t_new = new IndexedBlockDatatype(this->Count, this->Blocklen, &(this->Displ[0]), this->Basetype);

    return t_new;

}

void IndexedBlockDatatype::Codegen(Value *compactbuf, Value *scatteredbuf, Value* incount, bool pack) {

    Function* TheFunction = Builder.GetInsertBlock()->getParent();

    // Base address of the input buffer
    Value* scatteredbuf_orig_int = Builder.CreatePtrToInt(scatteredbuf, LLVM_INT64);
    Value* extend = constNode((long)this->getExtent());
    Value* incount_64 = Builder.CreateZExt(incount, LLVM_INT64);
    Value* incount_expanded = Builder.CreateMul(incount_64, extend);

    // Loop
    BasicBlock* PreheaderBB = Builder.GetInsertBlock();
    BasicBlock* LoopBB = BasicBlock::Create(getGlobalContext(), "loop", TheFunction);
    Builder.CreateBr(LoopBB);
    Builder.SetInsertPoint(LoopBB);

    PHINode *compact = Builder.CreatePHI(LLVM_INT8PTR, 2, "compact");
    compact->addIncoming(compactbuf, PreheaderBB);
    PHINode* i = Builder.CreatePHI(LLVM_INT64, 2, "i");
    i->addIncoming(constNode(0l), PreheaderBB);

    Value* compact_addr = Builder.CreatePtrToInt(compact, LLVM_INT64);

    // OPT: Make this the loop counter
    Value* scattered_disp_base = Builder.CreateAdd(scatteredbuf_orig_int, i);

    Value* nextcompact = compact;
    Value* compact_bytes_to_stride = constNode((long) Basetype->getSize() * Blocklen);

    for (int i=0; i<this->Count; i++) {
        // Set the scattered ptr to scattered_disp_base + this->Disl[i] * Basetype->size
        Value* displ_i = ConstantInt::get(getGlobalContext(), APInt(64, this->Displ[i] * Basetype->getSize(), false));
        Value* scattered_disp = Builder.CreateAdd(scattered_disp_base, displ_i);
        Value* scattered = Builder.CreateIntToPtr(scattered_disp, LLVM_INT8PTR);

        if (pack) Basetype->Codegen_Pack(scattered, ConstantInt::get(getGlobalContext(), APInt(32, Blocklen, false)), nextcompact);
        else      Basetype->Codegen_Unpack(nextcompact, ConstantInt::get(getGlobalContext(), APInt(32, Blocklen, false)), scattered);

        // Increment the compact ptr by Size(Basetype) * Blocklen
        compact_addr = Builder.CreateAdd(compact_addr, compact_bytes_to_stride);
        nextcompact = Builder.CreateIntToPtr(compact_addr, LLVM_INT8PTR);
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

void IndexedBlockDatatype::Codegen_Pack(Value* inbuf, Value* incount, Value* outbuf) {
    Codegen(outbuf, inbuf, incount, true);
}

void IndexedBlockDatatype::Codegen_Unpack(Value* inbuf, Value* incount, Value* outbuf) {
    Codegen(inbuf, outbuf, incount, false);
}

int IndexedBlockDatatype::getExtent() {

    if (this->Count == 0) return 0;

    int bext = this->Basetype->getExtent();

    int ub = this->Displ[0] * bext + this->Blocklen * bext;
    int lb = this->Displ[0] * bext;

    for (int i=0; i<this->Count; i++) {
        int tmp_ub = this->Displ[i] * bext + this->Blocklen * bext;
        int tmp_lb = this->Displ[i] * bext;
        if (tmp_ub > ub) ub = tmp_ub;
        if (tmp_lb < lb) lb = tmp_lb;
    }

    return ub - lb;

}

int IndexedBlockDatatype::getSize() {

    int sum = 0;
    int bsize = this->Basetype->getSize();
    for (int i=0; i<this->Count; i++) {
        sum += bsize * this->Blocklen;
    }

    return sum;

}

string IndexedBlockDatatype::toString() {
	stringstream res;
	res << "idxb(" << this->Count << " " << this->Blocklen << ")[";
    for (unsigned int i=0; i<Displ.size(); i++) {
		res << Displ[i];
		if (i <Displ.size() - 1) {
			res << " ";
		}
    }
	res << "]{" << Basetype->toString() << "}";
	return res.str();
}


/* Class HIndexedDatatype */
HIndexedDatatype::HIndexedDatatype(int count, int* blocklen, long* displ, Datatype* basetype) : Datatype() {

    this->Count = count;
    this->Basetype = basetype->Clone();
    for (int i=0; i<count; i++) this->Blocklen.push_back(blocklen[i]);
    for (int i=0; i<count; i++) this->Displ.push_back(displ[i]);

}

HIndexedDatatype::~HIndexedDatatype(void) {

    delete(this->Basetype);

}

HIndexedDatatype* HIndexedDatatype::Clone() {

    HIndexedDatatype* t_new = new HIndexedDatatype(this->Count, &(this->Blocklen[0]), &(this->Displ[0]), this->Basetype);

    return t_new;

}

void HIndexedDatatype::Codegen(Value *compactbuf, Value *scatteredbuf, Value* incount, bool pack) {
    Function* TheFunction = Builder.GetInsertBlock()->getParent();

    // Base address of the input buffer
    Value* scatteredbuf_orig_int = Builder.CreatePtrToInt(scatteredbuf, LLVM_INT64);
    Value* extend = constNode((long)this->getExtent());
    Value* incount_64 = Builder.CreateZExt(incount, LLVM_INT64);
    Value* incount_expanded = Builder.CreateMul(incount_64, extend);

    // Loop
    BasicBlock* PreheaderBB = Builder.GetInsertBlock();
    BasicBlock* LoopBB = BasicBlock::Create(getGlobalContext(), "loop", TheFunction);
    Builder.CreateBr(LoopBB);
    Builder.SetInsertPoint(LoopBB);

    PHINode *compact = Builder.CreatePHI(LLVM_INT8PTR, 2, "compact");
    compact->addIncoming(compactbuf, PreheaderBB);
    PHINode* i = Builder.CreatePHI(LLVM_INT64, 2, "i");
    i->addIncoming(constNode(0l), PreheaderBB);

    Value* compact_addr = Builder.CreatePtrToInt(compact, LLVM_INT64);

    // OPT: Make this the loop counter
    Value* scattered_disp_base = Builder.CreateAdd(scatteredbuf_orig_int, i);

    Value* nextcompact = compact;
    for (int i=0; i<this->Count; i++) {
        // Set the scattered ptr to scattered_disp_base + this->Disl[i]
        Value* displ_i = ConstantInt::get(getGlobalContext(), APInt(64, this->Displ[i], false));
        Value* scattered_disp = Builder.CreateAdd(scattered_disp_base, displ_i);
        Value* scattered = Builder.CreateIntToPtr(scattered_disp, LLVM_INT8PTR);

        if (pack) Basetype->Codegen_Pack(scattered, ConstantInt::get(getGlobalContext(), APInt(32, this->Blocklen[i], false)), nextcompact);
        else      Basetype->Codegen_Unpack(nextcompact, ConstantInt::get(getGlobalContext(), APInt(32, this->Blocklen[i], false)), scattered);

        // Increment the compact ptr by Size(Basetype) * Blocklen
        Value* compact_bytes_to_stride = constNode((long)Basetype->getSize() * this->Blocklen[i]);
        compact_addr = Builder.CreateAdd(compact_addr, compact_bytes_to_stride);
        nextcompact = Builder.CreateIntToPtr(compact_addr, LLVM_INT8PTR);
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

void HIndexedDatatype::Codegen_Pack(Value* inbuf, Value* incount, Value* outbuf) {
    Codegen(outbuf, inbuf, incount, true);
}

void HIndexedDatatype::Codegen_Unpack(Value* inbuf, Value* incount, Value* outbuf) {
    Codegen(inbuf, outbuf, incount, false);
}

int HIndexedDatatype::getExtent() {

    if (this->Count == 0) return 0;

    int bext = this->Basetype->getExtent();

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

int HIndexedDatatype::getSize() {

    int sum = 0;
    int bsize = this->Basetype->getSize();
    for (int i=0; i<this->Count; i++) {
        sum += bsize * this->Blocklen[i];
    }

    return sum;

}

string HIndexedDatatype::toString() {
	stringstream res;
	res << "hidx(";
    for (unsigned int i=0; i<Displ.size(); i++) {
		res << Displ[i] << "," << Blocklen[i];
		if (i <Displ.size() - 1) {
			res << " ";
		}
    }
	res << ")[" << Basetype->toString() << "]";
	return res.str();
}


/* Class StructDatatype */
StructDatatype::StructDatatype(int count, int* blocklen, long*  displ, Datatype** types) : Datatype() {

    this->Count = count;
    for (int i=0; i<count; i++) Blocklen.push_back(blocklen[i]);
    for (int i=0; i<count; i++) Displ.push_back(displ[i]);
    for (int i=0; i<count; i++) Types.push_back(types[i]->Clone());

}

StructDatatype::~StructDatatype(void) {

    for (int i=0; i<this->Count; i++) delete(Types[i]);

}

StructDatatype* StructDatatype::Clone() {

    StructDatatype* t_new = new StructDatatype(this->Count, &(this->Blocklen[0]), &(this->Displ[0]), &(this->Types[0]));

    return t_new;

}

void StructDatatype::Codegen(Value *compactbuf, Value *scatteredbuf, Value* incount, bool pack) {

    Function* TheFunction = Builder.GetInsertBlock()->getParent();

    // Base address of the input buffer
    Value* scatteredbuf_orig_int = Builder.CreatePtrToInt(scatteredbuf, LLVM_INT64);
    Value* extend = constNode((long)this->getExtent());
    Value* incount_64 = Builder.CreateZExt(incount, LLVM_INT64);
    Value* incount_expanded = Builder.CreateMul(incount_64, extend);

    // Loop
    BasicBlock* PreheaderBB = Builder.GetInsertBlock();
    BasicBlock* LoopBB = BasicBlock::Create(getGlobalContext(), "loop", TheFunction);
    Builder.CreateBr(LoopBB);
    Builder.SetInsertPoint(LoopBB);

    PHINode *compact = Builder.CreatePHI(LLVM_INT8PTR, 2, "compact");
    compact->addIncoming(compactbuf, PreheaderBB);
    PHINode* i = Builder.CreatePHI(LLVM_INT64, 2, "i");
    i->addIncoming(constNode(0l), PreheaderBB);

    Value* compact_addr = Builder.CreatePtrToInt(compact, LLVM_INT64);

    // OPT: Make this the loop counter
    Value* scattered_disp_base = Builder.CreateAdd(scatteredbuf_orig_int, i);

    Value* nextcompact = compact;
    for (int i=0; i<this->Count; i++) {
        // Set the scattered ptr to scattered_disp_base + this->Disl[i]
        Value* displ_i = ConstantInt::get(getGlobalContext(), APInt(64, this->Displ[i], false));
        Value* scattered_disp = Builder.CreateAdd(scattered_disp_base, displ_i);
        Value* scattered = Builder.CreateIntToPtr(scattered_disp, LLVM_INT8PTR);

        if (pack) Types[i]->Codegen_Pack(scattered, ConstantInt::get(getGlobalContext(), APInt(32, this->Blocklen[i], false)), nextcompact);
        else      Types[i]->Codegen_Unpack(nextcompact, ConstantInt::get(getGlobalContext(), APInt(32, this->Blocklen[i], false)), scattered);

        // Increment the compact ptr by Size(Basetype) * Blocklen
        Value* compact_bytes_to_stride = constNode((long)Types[i]->getSize() * this->Blocklen[i]);
        compact_addr = Builder.CreateAdd(compact_addr, compact_bytes_to_stride);
        nextcompact = Builder.CreateIntToPtr(compact_addr, LLVM_INT8PTR);
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

void StructDatatype::Codegen_Pack(Value* inbuf, Value* incount, Value* outbuf) {
    Codegen(outbuf, inbuf, incount, true);
}

void StructDatatype::Codegen_Unpack(Value* inbuf, Value* incount, Value* outbuf) {
    Codegen(inbuf, outbuf, incount, false);
}

int StructDatatype::getExtent() {

    if (this->Count == 0) return 0;

    int lb = this->Displ[0];
    int ub = this->Displ[0] + this->Types[0]->getExtent() * this->Blocklen[0];
    for (int i=0; i<this->Count; i++) {
        int tmp_ub = this->Displ[i] + this->Types[i]->getExtent() * this->Blocklen[i];
        int tmp_lb = this->Displ[i];
        if (tmp_ub > ub) ub = tmp_ub;
        if (tmp_lb < lb) lb = tmp_lb;
    }

    return ub-lb;

}

int StructDatatype::getSize() {

    int sum = 0;
    for (int i=0; i<this->Count; i++) {
        sum += this->Types[i]->getSize() * this->Blocklen[i];
    }

    return sum;

}

string StructDatatype::toString() {
	stringstream res;
	res << "struct(" << this->Count << ")";
    for (unsigned int i=0; i<Displ.size(); i++) {
		res << "[" << this->Displ[i] << "," << this->Blocklen[i] << "]";
		res << "{" << this->Types[i]->toString() << "}";
		if (i <Displ.size() - 1) {
			res << " ";
		}
    }
	return res.str();
}


// this jits the pack/unpack functions
void generate_pack_function(Datatype* ddt) {
#if TIME
    HRT_GET_TIMESTAMP(start);     
#endif

    Function* F = Function::Create(FT, Function::ExternalLinkage, "packer", TheModule);
    F->setDoesNotThrow();
    F->setDoesNotAlias(1);
    F->setDoesNotAlias(3);

    // Set names for all arguments.
    unsigned Idx = 0;
    for (Function::arg_iterator AI = F->arg_begin(); Idx != Args.size(); ++AI, ++Idx) {
		assert(AI !=  F->arg_end());
        AI->setName(Args[Idx]);
        NamedValues[Args[Idx]] = AI;
    }

    // Create a new basic block to start insertion into.
    BasicBlock *BB = BasicBlock::Create(getGlobalContext(), "entry", F);
    Builder.SetInsertPoint(BB);

    // generate code for the datatype
    ddt->Codegen_Pack(NamedValues["inbuf"], NamedValues["count"], NamedValues["outbuf"]);
    Builder.CreateRetVoid();

#if LLVM_VERIFY
    //F->viewCFG();
    verifyFunction(*F);
#endif
#if LLVM_OPTIMIZE
    TheFPM->run(*F);
#endif
#if LLVM_OUTPUT
    F->dump();

    std::vector<Type *> arg_type;
    arg_type.push_back(LLVM_INT8PTR);
    arg_type.push_back(LLVM_INT8PTR);
    arg_type.push_back(LLVM_INT64);
    Function *memcopy = Intrinsic::getDeclaration(TheModule, Intrinsic::memcpy, arg_type);
    memcopy->dump();

//    std::vector<Type *> prefetch_arg_type;
//    Function *prefetch = Intrinsic::getDeclaration(TheModule, Intrinsic::prefetch, prefetch_arg_type);
//    prefetch->dump();
#endif

    ddt->packer = (void (*)(void*, int, void*))(intptr_t) TheExecutionEngine->getPointerToFunction(F);
	ddt->FPack = F;

#if TIME
    HRT_GET_TIMESTAMP(stop);
    HRT_GET_ELAPSED_TICKS(start, stop, &tmp);
    commit_time += HRT_GET_USEC(tmp);
#endif
}

void generate_unpack_function(Datatype* ddt) {
#if TIME
    HRT_GET_TIMESTAMP(start);     
#endif

    Function* F = Function::Create(FT, Function::ExternalLinkage, "unpacker", TheModule);
    F->setDoesNotThrow();
    F->setDoesNotAlias(1);
    F->setDoesNotAlias(3);

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
    ddt->Codegen_Unpack(NamedValues["inbuf"], NamedValues["count"], NamedValues["outbuf"]);
    Builder.CreateRetVoid();

#if LLVM_VERIFY
    verifyFunction(*F);
#endif
#if LLVM_OPTIMIZE
    TheFPM->run(*F);
#endif
#if LLVM_OUTPUT
    F->dump();

    std::vector<Type *> arg_type;
    arg_type.push_back(LLVM_INT8PTR);
    arg_type.push_back(LLVM_INT8PTR);
    arg_type.push_back(LLVM_INT64);
    Function *memcopy = Intrinsic::getDeclaration(TheModule, Intrinsic::memcpy, arg_type);
    memcopy->dump();

//    std::vector<Type *> prefetch_arg_type;
//    Function *prefetch = Intrinsic::getDeclaration(TheModule, Intrinsic::prefetch, prefetch_arg_type);
//    prefetch->dump();
#endif

    ddt->unpacker = (void (*)(void*, int, void*))(intptr_t) TheExecutionEngine->getPointerToFunction(F);
    ddt->FUnpack = F;

#if TIME
    HRT_GET_TIMESTAMP(stop);
    HRT_GET_ELAPSED_TICKS(start, stop, &tmp);
    commit_time += HRT_GET_USEC(tmp);
#endif
}

void DDT_Commit(Datatype* ddt) {
#if DDT_OUTPUT
    ddt->print("");
#endif
#if !LAZY
    generate_pack_function(ddt);
    generate_unpack_function(ddt);
#endif
}

// this calls the pack/unpack function
void DDT_Pack(void* inbuf, void* outbuf, Datatype* ddt, int count) {
#if LAZY
    if (ddt->packer == NULL) generate_pack_function(ddt);
#endif
    ddt->packer(inbuf, count, outbuf);
}

void DDT_Lazy_Unpack_Commit(Datatype* ddt) {
#if LAZY
    if (ddt->unpacker == NULL) generate_unpack_function(ddt);
#endif
}

void DDT_Unpack(void* inbuf, void* outbuf, Datatype* ddt, int count) {
#if LAZY
    DDT_Lazy_Unpack_Commit(ddt);
#endif
    ddt->unpacker(inbuf, count, outbuf);
}

void DDT_Free(Datatype* ddt) {

	if (ddt->packer != NULL) {
		TheExecutionEngine->freeMachineCodeForFunction(ddt->FPack);
		ddt->FPack->eraseFromParent();
	}
	if (ddt->unpacker != NULL) {
		TheExecutionEngine->freeMachineCodeForFunction(ddt->FUnpack);
		ddt->FUnpack->eraseFromParent();
	}
    delete ddt;
    
}

// init the JIT compiler
void DDT_Init() {
#if TIME
    HRT_INIT(1, g_timerfreq);
#endif

    InitializeNativeTarget();
    LLVMContext &Context = getGlobalContext();
    TheModule = new Module("FARC-JIT", Context);

    // Create the JIT.  This takes ownership of the module.
    std::string ErrStr;
    EngineBuilder engine_builder(TheModule);
    engine_builder.setEngineKind(EngineKind::JIT);
    engine_builder.setOptLevel(CodeGenOpt::Aggressive);
    engine_builder.setErrorStr(&ErrStr);

    TheExecutionEngine = engine_builder.create();

    if (!TheExecutionEngine) {
        fprintf(stderr, "Could not create ExecutionEngine: %s\n", ErrStr.c_str());
        exit(1);
    }

    // Initialize some types used by all packers
    std::vector<Type*> FuncArgs;
    FuncArgs.push_back(LLVM_INT8PTR);
    FuncArgs.push_back(LLVM_INT32);
    FuncArgs.push_back(LLVM_INT8PTR);
    FT = FunctionType::get(LLVM_VOID, FuncArgs, false);

    Args.push_back("inbuf");
    Args.push_back("count");
    Args.push_back("outbuf");


#if LLVM_OPTIMIZE
    FunctionPassManager* OurFPM = new FunctionPassManager(TheModule);

    /*
    PassManagerBuilder Builder;
    Builder.OptLevel = 3;
    Builder.Vectorize = true;
    Builder.LoopVectorize = true;
    Builder.populateFunctionPassManager(*OurFPM);
    */

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
#endif

}

void DDT_Finalize() {
#if TIME
    printf("Commit time: %10.3lf s\n", commit_time/1000000);
#endif
}

} // namespace farc

/* Define the functions of the pack header file */
