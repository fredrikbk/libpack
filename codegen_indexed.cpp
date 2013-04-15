#include "codegen.hpp"
#include "codegen_common.hpp"
#include "ddt_jit.hpp"

#include <llvm/IR/Value.h>

using namespace llvm;
using namespace std;

namespace farc {

#define IDXB_LOOP_TRESHOLD  16
#define IDXB_LOOP_UNROLL    1

void codegenIndexedBlock(Value *compactbuf, Value *scatteredbuf, Value* incount,
                         int extent, int size, int count, int blocklen, Datatype *basetype,
                         const vector<int> &displs, Value* indices_arr,
                         bool pack) {
    Function* func = Builder.GetInsertBlock()->getParent();

    if (count > IDXB_LOOP_TRESHOLD) {
        // Entry block
        BasicBlock* preamble = Builder.GetInsertBlock();
        Value* noncontig = Builder.CreatePtrToInt(scatteredbuf, LLVM_INT64);
        noncontig->setName("noncontig");
        Value* contig = Builder.CreatePtrToInt(compactbuf, LLVM_INT64);
        contig->setName("contig");

        Value* incount64   = Builder.CreateZExt(incount, LLVM_INT64, "count64");
        Value* bytestocopy = Builder.CreateMul(incount64, constNode((long)size), "bytestocopy");
        Value* exitcond = Builder.CreateAdd(contig, bytestocopy, "exitcond");

        // Outer loop
        BasicBlock* outerloop = BasicBlock::Create(getGlobalContext(), "outerloop", func);
        Builder.CreateBr(outerloop);
        Builder.SetInsertPoint(outerloop);
        
        PHINode* contig1 = Builder.CreatePHI(LLVM_INT64, 2, "contig1");
        contig1->addIncoming(contig, preamble);
        PHINode* noncontig1 = Builder.CreatePHI(LLVM_INT64, 2, "noncontig1");
        noncontig1->addIncoming(noncontig, preamble);

        // Inner loop
        BasicBlock *innerloop = BasicBlock::Create(getGlobalContext(), "innerloop", func);
        Builder.CreateBr(innerloop);
        Builder.SetInsertPoint(innerloop);

        PHINode *i = Builder.CreatePHI(LLVM_INT64, 2, "i");
        i->addIncoming(constNode((long)0), outerloop);
        PHINode* contig2 = Builder.CreatePHI(LLVM_INT64, 2, "contig2");
        contig2->addIncoming(contig1, outerloop);                

        std::vector<Value*> arrayidx_list;
        arrayidx_list.push_back(constNode((long)0));
        arrayidx_list.push_back(i);
        Value* displloc = Builder.CreateGEP(indices_arr, arrayidx_list, "displloc");
        Value* displ = Builder.CreateLoad(displloc, "displ");
        Value* displ64 = Builder.CreateZExt(displ, LLVM_INT64, "displ64");
        Value *noncontig2 = Builder.CreateAdd(noncontig1, displ64, "noncontig2");

        Value* contig2ptr = Builder.CreateIntToPtr(contig2, LLVM_INT8PTR, "contig2ptr");
        Value* noncontig2ptr = Builder.CreateIntToPtr(noncontig2, LLVM_INT8PTR, "noncontig2ptr");

        if (pack) basetype->packCodegen(noncontig2ptr, constNode(blocklen), contig2ptr);
        else      basetype->unpackCodegen(contig2ptr, constNode(blocklen), noncontig2ptr);

        Value* nextcontig2 =
            Builder.CreateAdd(contig2, constNode((long)blocklen*basetype->getExtent()), "nextcontig2");
        Value* nexti = Builder.CreateAdd(i, constNode((long)IDXB_LOOP_UNROLL), "nexti");

        contig2->addIncoming(nextcontig2, innerloop);
        i->addIncoming(nexti, innerloop);
                
        Value* innertest =
            Builder.CreateICmpEQ(nexti, constNode((long)count), "innertest");
        BasicBlock *innerpost = BasicBlock::Create(getGlobalContext(), "innerpost", func);
        Builder.CreateCondBr(innertest, innerpost, innerloop);
        Builder.SetInsertPoint(innerpost);
        // End of inner loop

        Value* nextnoncontig1 =
            Builder.CreateAdd(noncontig1, constNode((long)extent), "nextnoncontig1");
        Value* nextcontig1 =
            Builder.CreateAdd(contig1, constNode((long)size), "nextcontig1");

        noncontig1->addIncoming(nextnoncontig1, innerpost);
        contig1->addIncoming(nextcontig1, innerpost);

        Value* outertest = Builder.CreateICmpEQ(nextcontig1, exitcond, "outertest");
        BasicBlock *outerpost = BasicBlock::Create(getGlobalContext(), "outerpost", func);
        Builder.CreateCondBr(outertest, outerpost, outerloop);
        Builder.SetInsertPoint(outerpost);
        // End of outer loop
    }
    else {
        // Base address of the input buffer
        Value* scatteredbuf_orig_int = Builder.CreatePtrToInt(scatteredbuf, LLVM_INT64);
        Value* extend = constNode((long)extent);
        Value* incount_64 = Builder.CreateZExt(incount, LLVM_INT64);
        Value* incount_expanded = Builder.CreateMul(incount_64, extend);

        // Loop
        BasicBlock* PreheaderBB = Builder.GetInsertBlock();
        BasicBlock* LoopBB = BasicBlock::Create(getGlobalContext(), "loop", func);
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
        Value* compact_bytes_to_stride = constNode((long)basetype->getSize() * blocklen);

        for (int i=0; i<count; i++) {
            // Set the scattered ptr to scattered_disp_base + this->Disl[i] * Basetype->size
            Value* displ_i = constNode((long)displs[i] * basetype->getSize());
            Value* scattered_disp = Builder.CreateAdd(scattered_disp_base, displ_i);
            Value* scattered = Builder.CreateIntToPtr(scattered_disp, LLVM_INT8PTR);

            if (pack) basetype->packCodegen(scattered, constNode(blocklen), nextcompact);
            else      basetype->unpackCodegen(nextcompact, constNode(blocklen), scattered);

            // Increment the compact ptr by Size(Basetype) * Blocklen
            compact_addr = Builder.CreateAdd(compact_addr, compact_bytes_to_stride);
            nextcompact = Builder.CreateIntToPtr(compact_addr, LLVM_INT8PTR);
        }

        // Increment the loop index and test for loop exit
        Value* nexti = Builder.CreateAdd(i, extend, "nexti");
        Value* EndCond = Builder.CreateICmpEQ(nexti, incount_expanded, "loopcond");

        // Create and branch to the outer loop postamble
        BasicBlock *LoopEndBB = Builder.GetInsertBlock();
        BasicBlock *AfterBB = BasicBlock::Create(getGlobalContext(), "afterloop", func);

        Builder.CreateCondBr(EndCond, AfterBB, LoopBB);
        Builder.SetInsertPoint(AfterBB);
                            
        // Add backedges for the loop induction variable
        compact->addIncoming(nextcompact, LoopEndBB);
        i->addIncoming(nexti, LoopEndBB);
    }
}

void codegenHindexed(Value *compactbuf, Value *scatteredbuf, Value* incount,
                     int extent, int count, Datatype *basetype,
                     const vector<int> &blocklens, const vector<long> &displs,
                     bool pack) {

    Function* func = Builder.GetInsertBlock()->getParent();

    // Base address of the input buffer
    Value* scatteredbuf_orig_int = Builder.CreatePtrToInt(scatteredbuf, LLVM_INT64);
    Value* extend = constNode((long)extent);
    Value* incount_64 = Builder.CreateZExt(incount, LLVM_INT64);
    Value* incount_expanded = Builder.CreateMul(incount_64, extend);

    // Loop
    BasicBlock* PreheaderBB = Builder.GetInsertBlock();
    BasicBlock* LoopBB = BasicBlock::Create(getGlobalContext(), "loop", func);
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
    for (int i=0; i<count; i++) {
        // Set the scattered ptr to scattered_disp_base + this->Disl[i]
        Value* displ_i = constNode((long)displs[i]);
        Value* scattered_disp = Builder.CreateAdd(scattered_disp_base, displ_i);
        Value* scattered = Builder.CreateIntToPtr(scattered_disp, LLVM_INT8PTR);

        if (pack) basetype->packCodegen(scattered, constNode(blocklens[i]), nextcompact);
        else      basetype->unpackCodegen(nextcompact, constNode(blocklens[i]), scattered);

        // Increment the compact ptr by Size(Basetype) * Blocklen
        Value* compact_bytes_to_stride = constNode((long)basetype->getSize() * blocklens[i]);
        compact_addr = Builder.CreateAdd(compact_addr, compact_bytes_to_stride);
        nextcompact = Builder.CreateIntToPtr(compact_addr, LLVM_INT8PTR);
    }

    // Increment the loop index and test for loop exit
    Value* nexti = Builder.CreateAdd(i, extend, "nexti");
    Value* EndCond = Builder.CreateICmpEQ(nexti, incount_expanded, "loopcond");

    // Create and branch to the outer loop postamble
    BasicBlock *LoopEndBB = Builder.GetInsertBlock();
    BasicBlock *AfterBB = BasicBlock::Create(getGlobalContext(), "afterloop", func);

    Builder.CreateCondBr(EndCond, AfterBB, LoopBB);
    Builder.SetInsertPoint(AfterBB);

    // Add backedges for the loop induction variable
    compact->addIncoming(nextcompact, LoopEndBB);
    i->addIncoming(nexti, LoopEndBB);
}

void codegenStruct(Value *compactbuf, Value *scatteredbuf,
                   Value* incount, int extent, int count,
                   const vector<int> &blocklens,
                   const vector<long> &displs,
                   const vector<Datatype*> &basetypes,
                   bool pack) {

    Function* func = Builder.GetInsertBlock()->getParent();

    // Base address of the input buffer
    Value* scatteredbuf_orig_int = Builder.CreatePtrToInt(scatteredbuf, LLVM_INT64);
    Value* extend = constNode((long)extent);
    Value* incount_64 = Builder.CreateZExt(incount, LLVM_INT64);
    Value* incount_expanded = Builder.CreateMul(incount_64, extend);

    // Loop
    BasicBlock* PreheaderBB = Builder.GetInsertBlock();
    BasicBlock* LoopBB = BasicBlock::Create(getGlobalContext(), "loop", func);
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
    for (int i=0; i<count; i++) {
        // Set the scattered ptr to scattered_disp_base + this->Disl[i]
        Value* displ_i = ConstantInt::get(getGlobalContext(), APInt(64, displs[i], false));
        Value* scattered_disp = Builder.CreateAdd(scattered_disp_base, displ_i);
        Value* scattered = Builder.CreateIntToPtr(scattered_disp, LLVM_INT8PTR);

        if (pack) basetypes[i]->packCodegen(scattered, constNode(blocklens[i]), nextcompact);
        else      basetypes[i]->unpackCodegen(nextcompact, constNode(blocklens[i]), scattered);

        // Increment the compact ptr by Size(Basetype) * Blocklen
        Value* compact_bytes_to_stride =
            constNode((long)basetypes[i]->getSize() * blocklens[i]);
        compact_addr = Builder.CreateAdd(compact_addr, compact_bytes_to_stride);
        nextcompact = Builder.CreateIntToPtr(compact_addr, LLVM_INT8PTR);
    }

    // Increment the loop index and test for loop exit
    Value* nexti = Builder.CreateAdd(i, extend, "nexti");
    Value* EndCond = Builder.CreateICmpEQ(nexti, incount_expanded, "loopcond");

    // Create and branch to the outer loop postamble
    BasicBlock *LoopEndBB = Builder.GetInsertBlock();
    BasicBlock *AfterBB = BasicBlock::Create(getGlobalContext(), "afterloop", func);

    Builder.CreateCondBr(EndCond, AfterBB, LoopBB);
    Builder.SetInsertPoint(AfterBB);
                            
    // Add backedges for the loop induction variable
    compact->addIncoming(nextcompact, LoopEndBB);
    i->addIncoming(nexti, LoopEndBB);
}

}
