// Copyright 2013 Timo Schneider and Fredrik Berg Kjolstad
//
// This file is part of the libpack packing library.
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT in the top level directory for details.

#ifndef CODEGEN_COMMON_H
#define CODEGEN_COMMON_H

#include "ddt_jit.hpp"

#include <llvm/IR/Value.h>
#include <llvm/IR/Constants.h>
#include <llvm/IR/IRBuilder.h>

#define LLVM_VOID     Type::getVoidTy(getGlobalContext())
#define LLVM_INT      Type::getInt32Ty(getGlobalContext())
#define LLVM_INT8     Type::getInt8Ty(getGlobalContext())
#define LLVM_INT32    Type::getInt32Ty(getGlobalContext())
#define LLVM_INT64    Type::getInt64Ty(getGlobalContext())
#define LLVM_INT8PTR  Type::getInt8PtrTy(getGlobalContext())
#define LLVM_FLOAT    Type::getFloatTy(getGlobalContext())
#define LLVM_DOUBLE   Type::getDoubleTy(getGlobalContext())

namespace farc {

extern llvm::IRBuilder<> Builder;

llvm::Value* multNode(int op1, llvm::Value* op2PtrNode);
llvm::ConstantInt* constNode(int val);
llvm::ConstantInt* constNode(long val);

void vmove(llvm::Value *dst, llvm::Value *src, int count, llvm::Type *elemtype);
llvm::Value *incrementPtr(llvm::Value *ptr, int byteInc);
llvm::Type *toLLVMType(PrimitiveDatatype::PrimitiveType type);

}

#endif // CODEGEN_COMMON_H
