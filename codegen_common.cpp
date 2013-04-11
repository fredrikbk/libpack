#include "codegen_common.hpp"

#include <cstdio>
#include <llvm/IR/Constants.h>
#include <llvm/IR/LLVMContext.h>

using namespace llvm;

namespace farc {

IRBuilder<> Builder(getGlobalContext());

Value* multNode(int op1, Value* op2PtrNode) {
    Value* op1Node = constNode((long)op1);
    Value* op2Node = Builder.CreateIntCast(op2PtrNode, LLVM_INT64, false); 
    return Builder.CreateMul(op1Node, op2Node);
}

Value* constNode(int val) {
    return ConstantInt::get(getGlobalContext(), APInt(32, val, false));
}

Value* constNode(long val) {
    return ConstantInt::get(getGlobalContext(), APInt(64, val, false));
}

void vmove(Value *dst, Value *src, int count, Type *elemtype) {
	Type *elemvectype_ptr = PointerType::getUnqual(VectorType::get(elemtype, count));
	Value *in_vec = Builder.CreateBitCast(src, elemvectype_ptr, "in2_addr_vec");
	Value *out_vec = Builder.CreateBitCast(dst, elemvectype_ptr, "out2_addr_vec");
	Value *elems = Builder.CreateAlignedLoad(in_vec, 1, "elems");
	Builder.CreateAlignedStore(elems, out_vec, 1);
}

Value *incrementPtr(Value *ptr, int byteInc) {
	Value *addr = Builder.CreatePtrToInt(ptr, LLVM_INT64);
	Value *newaddr = Builder.CreateAdd(addr, Builder.getInt64(byteInc));
	return Builder.CreateIntToPtr(newaddr, LLVM_INT8PTR);
}

Type *toLLVMType(PrimitiveDatatype::PrimitiveType type) {
	Type *elemtype = NULL;
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

}
