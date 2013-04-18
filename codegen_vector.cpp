// Copyright 2013 Timo Schneider and Fredrik Berg Kjolstad
//
// This file is part of the libpack packing library.
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT in the top level directory for details.

#include "codegen.hpp"
#include "codegen_common.hpp"
#include "ddt_jit.hpp"

#include <llvm/IR/Value.h>

using namespace llvm;

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

namespace farc {

void codegenVector(Value* inbuf, Value* incount, Value* outbuf,
                   Datatype* basetype, int count, int blocklen,
                   int elemstride_in, int elemstride_out, bool pack) {
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
    if (pack) basetype->packCodegen(in2_addr, ConstantInt::get(getGlobalContext(), APInt(32, blocklen, false)), out2_addr);
    else      basetype->unpackCodegen(in2_addr, ConstantInt::get(getGlobalContext(), APInt(32, blocklen, false)), out2_addr);


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

}
