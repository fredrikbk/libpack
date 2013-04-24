// Copyright 2013 Timo Schneider and Fredrik Berg Kjolstad
//
// This file is part of the libpack packing library.
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT in the top level directory for details.

#include "codegen.hpp"
#include "codegen_common.hpp"
#include "ddt_jit.hpp"
#include <stdio.h>

#include <llvm/IR/Value.h>

using namespace llvm;

namespace farc {

void codegenContiguous(Value* inbuf, Value* incount,
                       Value* outbuf, Datatype *basetype,
                       int inptr_inc, int outptr_inc,
                       int count, bool pack) {
    Function* TheFunction = Builder.GetInsertBlock()->getParent();

    // Loop
    BasicBlock* PreheaderBB = Builder.GetInsertBlock();
    BasicBlock* LoopBB =
        BasicBlock::Create(getGlobalContext(), "loop", TheFunction);
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
    if (pack) basetype->packCodegen(in, constNode(count), out);
    else      basetype->unpackCodegen(in, constNode(count), out);


    // Increment the out ptr by outptr_inc
    Value* out_bytes_to_stride = constNode((long) outptr_inc);
    Value* out_addr_cvi = Builder.CreatePtrToInt(out, LLVM_INT64);
    Value* out_addr = Builder.CreateAdd(out_addr_cvi, out_bytes_to_stride);
    Value* nextout = Builder.CreateIntToPtr(out_addr, LLVM_INT8PTR);

    // Increment the in ptr by inptr_inc
    Value* in_bytes_to_stride = constNode((long) inptr_inc);
    Value* in_addr_cvi = Builder.CreatePtrToInt(in, LLVM_INT64);
    Value* in_addr = Builder.CreateAdd(in_addr_cvi, in_bytes_to_stride);
    Value* nextin = Builder.CreateIntToPtr(in_addr, LLVM_INT8PTR);

    // Increment outer loop index
    Value* nexti = Builder.CreateAdd(i, constNode(1), "nexti");
    Value* EndCond_outer = Builder.CreateICmpEQ(nexti, incount, "loopcond");

    // Create and branch to the outer loop postamble
    BasicBlock *LoopEndBB = Builder.GetInsertBlock();
    BasicBlock *AfterBB =
        BasicBlock::Create(getGlobalContext(), "afterloop", TheFunction);
    Builder.CreateCondBr(EndCond_outer, AfterBB, LoopBB);
    Builder.SetInsertPoint(AfterBB);

    // Add backedges for the outer loop induction variable
    out->addIncoming(nextout, LoopEndBB);
    in->addIncoming(nextin, LoopEndBB);
    i->addIncoming(nexti, LoopEndBB);
}

}
