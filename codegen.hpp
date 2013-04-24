// Copyright 2013 Timo Schneider and Fredrik Berg Kjolstad
//
// This file is part of the libpack packing library.
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT in the top level directory for details.

#ifndef CODEGEN_H
#define CODEGEN_H

#include "ddt_jit.hpp"
#include <vector>

namespace llvm {
class Value;
class GlobalVariable;
}

#define IDXB_LOOP_TRESHOLD  16
#define IDXB_LOOP_UNROLL    1

namespace farc {

void codegenPrimitive(llvm::Value* inbuf, llvm::Value* incount,
                      llvm::Value* outbuf, int size, 
                      PrimitiveDatatype::PrimitiveType type);

void codegenPrimitiveResized(llvm::Value* inbuf, llvm::Value* incount,
                      llvm::Value* outbuf, int size, int extent, 
                      PrimitiveDatatype::PrimitiveType type);

void codegenContiguous(llvm::Value* inbuf, llvm::Value* incount,
                       llvm::Value* outbuf, Datatype *basetype,
                       int elemstride_in, int elemstride_out,
                       int count, bool pack);

void codegenVector(llvm::Value *inbuf, llvm::Value *incount,
                   llvm::Value *outbuf, Datatype *basetype, int count,
                   int blocklen, int elemstride_in, int elemstride_out,
                   int inptr_inc, int outptr_inc, bool pack);

void codegenIndexedBlock(llvm::Value *compactbuf, llvm::Value *scatteredbuf,
                         llvm::Value* incount, int extent, int size,
                         int count, int blocklen, Datatype *basetype,
                         const std::vector<int> &displs,
                         llvm::Value* indices_arr,
                         bool pack);

void codegenHindexed(llvm::Value *compactbuf, llvm::Value *scatteredbuf,
                     llvm::Value* incount, int extent, int count,
                     Datatype *basetype, const std::vector<int> &blocklens,
                     const std::vector<long> &displs, bool pack);

void codegenStruct(llvm::Value *compactbuf, llvm::Value *scatteredbuf,
                   llvm::Value* incount, int extent, int count,
                   const std::vector<int> &blocklens,
                   const std::vector<long> &displs,
                   const std::vector<Datatype*> &basetypes,
                   bool pack);

}

#endif // CODEGEN_H
