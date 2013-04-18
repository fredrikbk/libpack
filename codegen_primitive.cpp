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

#if !PACKVAR
#undef PACKVAR
#define PACKVAR 1
#endif

namespace farc {

void codegenPrimitive(Value* inbuf, Value* incount, Value* outbuf,
                      int size, PrimitiveDatatype::PrimitiveType type) {
    Function* TheFunction = Builder.GetInsertBlock()->getParent();
    llvm::ConstantInt* incount_ci = dyn_cast<llvm::ConstantInt>(incount);

    if (incount_ci == NULL) {
        Value* contig_extend = multNode(size, incount);
        Builder.CreateMemCpy(outbuf, inbuf, contig_extend, 1);
    }
    else {
            
#if PACKVAR == 1 // Unaligned loads and stores
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
        #define COPY_LOOP_TRESHOLD 64*8
        #endif

        // Number of times to unroll the copy loop
        #ifndef COPY_LOOP_UNROLL
        #define COPY_LOOP_UNROLL 16
        #endif

        const int LOOP_ELEM_TRESHOLD = COPY_LOOP_TRESHOLD / size;
        const int vector_size = SIMD_BYTE_SIZE / size;

        // Assert power-of-two
        assert((SIMD_BYTE_SIZE & (SIMD_BYTE_SIZE-1)) == 0);
        assert((vector_size & (vector_size-1)) == 0);

        assert((SIMD_BYTE_SIZE % size) == 0);

        llvm::Type *elemtype = toLLVMType(type);

        // Number of elements to copy
        int incount_val = incount_ci->getSExtValue();

        // Trunk to a multiple of the unroll factor
        const int vectors_to_copy =
            ((incount_val / vector_size) / COPY_LOOP_UNROLL) * COPY_LOOP_UNROLL;

        // Copy vectors_to_copy vectors of size vecsize.
        // If the number of values is less than the LOOP_ELEM_TRESHOLD
        // treshold then we fall through to the postamble and unroll
        // everything
        if (vectors_to_copy > 0 && incount_val >= LOOP_ELEM_TRESHOLD) {
            Value *inbuf_int = Builder.CreatePtrToInt(inbuf, LLVM_INT64);
            Value *exitval   = Builder.CreateAdd(inbuf_int, constNode((long)vectors_to_copy * SIMD_BYTE_SIZE), "exitval");

            BasicBlock *header = Builder.GetInsertBlock();
            BasicBlock *copyloop   =
                BasicBlock::Create(getGlobalContext(), "copyloop", TheFunction);
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
                in_addr = Builder.CreateAdd(in_addr_cvi, Builder.getInt64(size * vector_size));
                inbuf = Builder.CreateIntToPtr(in_addr, LLVM_INT8PTR);

                outbuf = incrementPtr(outbuf, size * vector_size);
            }
            incount_val -= vectors_to_copy * vector_size;

            inphi->addIncoming(inbuf, copyloop);
            outphi->addIncoming(outbuf, copyloop);

            // Create and jump to postamble
            BasicBlock *copypostamble =
                BasicBlock::Create(getGlobalContext(), "copypostamble", TheFunction);
            Value *exitcond = Builder.CreateICmpEQ(in_addr, exitval);
            Builder.CreateCondBr(exitcond, copypostamble, copyloop);
            Builder.SetInsertPoint(copypostamble);
        }

        // Copy postamble: copy the overflow elements that did not fit in full vector
        for (int vecsize=vector_size; vecsize > 0; vecsize /= 2) {
            const int veccount = incount_val / vecsize;
            for (int i=0; i<veccount; i++) {
                vmove(outbuf, inbuf, vecsize, elemtype);
                inbuf  = incrementPtr(inbuf, size * vecsize);
                outbuf = incrementPtr(outbuf, size * vecsize);
            }
            incount_val -= veccount * vecsize;
        }
#elif PACKVAR == 2 // memcopy
        Value* contig_extend = multNode(this->getSize(), incount);
        Value* memcopy = Builder.CreateMemCpy(outbuf, inbuf, contig_extend, 1);
#elif PACKVAR == 3// aligned loads and stores
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
#else
#error NO PACKVAR DEFINED
#endif

    }
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

}
