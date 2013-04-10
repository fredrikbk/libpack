#ifndef CODEGEN_H
#define CODEGEN_H

#include "ddt_jit.hpp"
#include <llvm/IR/Value.h>
#include <vector>

namespace farc {

void codegenPrimitive(llvm::Value* inbuf, llvm::Value* incount,
                      llvm::Value* outbuf, int size, 
                      PrimitiveDatatype::PrimitiveType type);

void codegenContiguous(llvm::Value* inbuf, llvm::Value* incount,
                       llvm::Value* outbuf, Datatype *basetype,
                       int elemstride_in, int elemstride_out,
                       int count, bool pack);

void codegenVector(llvm::Value *inbuf, llvm::Value *incount,
                   llvm::Value *outbuf, Datatype *basetype, int count,
                   int blocklen, int elemstride_in, int elemstride_out,
                   bool pack);

void codegenIndexedBlock(llvm::Value *compactbuf, llvm::Value *scatteredbuf,
                         llvm::Value* incount, int extent, int count,
                         int blocklen, Datatype *basetype,
                         const std::vector<int> &displs, bool pack);

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