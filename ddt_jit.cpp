#include "ddt_jit.hpp"
#include "pack.h"

#include "codegen.hpp"
#include "codegen_common.hpp"

#include <map>
#include <cstdio>
#include <iostream>
#include <sstream>

#include "llvm/IR/Module.h"
#include "llvm/IR/DerivedTypes.h"
#include "llvm/IR/Intrinsics.h"

#include "llvm/Support/TargetSelect.h"
#include "llvm/ExecutionEngine/ExecutionEngine.h"
#include "llvm/ExecutionEngine/JIT.h"


#define LAZY           0

#define DDT_OPTIMIZE   1
#define DDT_OUTPUT     0 

#define LLVM_OPTIMIZE  0

// LLVM_OUTPUT should be picked up from the environment by the build system
#ifndef LLVM_OUTPUT
#define LLVM_OUTPUT 0
#endif

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

using namespace std;
using namespace llvm;

namespace farc {

static Module *module;
static std::map<std::string, Value*> NamedValues;
static ExecutionEngine *TheExecutionEngine;

std::vector<std::string> Args;
FunctionType *FT;


/* Datatype */
Datatype::~Datatype() {
    if (this->pack != NULL) {
        TheExecutionEngine->freeMachineCodeForFunction(this->fpack);
        this->fpack->eraseFromParent();
        this->pack = NULL;
    }
    if (this->unpack != NULL) {
        TheExecutionEngine->freeMachineCodeForFunction(this->funpack);
        this->funpack->eraseFromParent();
        this->unpack = NULL;
    }
}


static inline Function* createFunctionHeader(const char *name) {
    Function* F = Function::Create(FT, Function::ExternalLinkage, name, module);
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

    return F;
}

static inline void postProcessFunction(Function *F) {
#if LLVM_VERIFY
    //F->viewCFG();
    verifyFunction(*F);
#endif
#if LLVM_OPTIMIZE
    TheFPM->run(*F);
#endif
}

void Datatype::compile(CompilationType type) {
    // Compress the datatype, by substituting datatypes for
    // equivalent, but more compact, datatypes
    #if DDT_OPTIMIZE
    Datatype *ddt = this->compress();
    #else
    Datatype *ddt = this;
    #endif

    // Create the global arrays needed by the datatypes
    ddt->globalCodegen(module);

    bool pack   = (type == PACK_UNPACK || type == PACK)   ? true : false;
    bool unpack = (type == PACK_UNPACK || type == UNPACK) ? true : false;
    if (pack) {
        this->fpack = createFunctionHeader("pack");

        // Create a new basic block to start insertion into.
        BasicBlock *BB = BasicBlock::Create(getGlobalContext(), "entry", this->fpack);
        Builder.SetInsertPoint(BB);

        // generate code for the datatype
        ddt->packCodegen(NamedValues["inbuf"], NamedValues["count"], NamedValues["outbuf"]);
        Builder.CreateRetVoid();

        postProcessFunction(this->fpack);

        #if !LLVM_OUTPUT
        this->pack = (void (*)(void*,int,void*))(intptr_t)
            TheExecutionEngine->getPointerToFunction(this->fpack);
        #endif
    }

    if (unpack) {
        this->funpack = createFunctionHeader("unpack");

        // Create a new basic block to start insertion into.
        BasicBlock *BB = BasicBlock::Create(getGlobalContext(), "entry", this->funpack);
        Builder.SetInsertPoint(BB);

        // generate code for the datatype
        ddt->unpackCodegen(NamedValues["inbuf"], NamedValues["count"], NamedValues["outbuf"]);
        Builder.CreateRetVoid();

        postProcessFunction(this->funpack);

        #if !LLVM_OUTPUT
        this->unpack = (void (*)(void*,int,void*))(intptr_t)
            TheExecutionEngine->getPointerToFunction(this->funpack);
        #endif
    }

    #if LLVM_OUTPUT
    // std::vector<Type *> arg_type;
    // arg_type.push_back(LLVM_INT8PTR);
    // arg_type.push_back(LLVM_INT8PTR);
    // arg_type.push_back(LLVM_INT64);
    // Function *memcopy = Intrinsic::getDeclaration(module, Intrinsic::memcpy, arg_type);
    // memcopy->dump();

    // std::vector<Type *> prefetch_arg_type;
    // Function *prefetch = Intrinsic::getDeclaration(module,Intrinsic::prefetch, prefetch_arg_type);
    // prefetch->dump();

    module->dump();

    if (pack) {
        this->pack = (void (*)(void*,int,void*))(intptr_t)
            TheExecutionEngine->getPointerToFunction(this->fpack);
    }
    if (unpack) {
        this->unpack = (void (*)(void*,int,void*))(intptr_t)
            TheExecutionEngine->getPointerToFunction(this->funpack);
    }
    #endif

    #if DDT_OPTIMIZE
    delete ddt;
    #endif
}

void Datatype::print(bool summary) {
    printf("%s\n", this->toString(summary).c_str());
}


/* PrimitiveDatatype */
PrimitiveDatatype::PrimitiveDatatype(PrimitiveDatatype::PrimitiveType type) : Datatype() {
    this->type = type;

    if (type == BYTE)   this->extent = 1;
    if (type == CHAR)   this->extent = 1;
    if (type == DOUBLE) this->extent = sizeof(double);
    if (type == FLOAT)  this->extent = sizeof(float);
    if (type == INT)    this->extent = sizeof(int);
    //TODO add more. Remember to also add them to the print function.

    this->size = this->extent;
}

PrimitiveDatatype* PrimitiveDatatype::clone() {
    return new PrimitiveDatatype(this->type);
}

void PrimitiveDatatype::packCodegen(Value* inbuf, Value* incount,
                                     Value* outbuf) {
    codegenPrimitive(inbuf, incount, outbuf, this->size, this->type);
}

void PrimitiveDatatype::unpackCodegen(Value* inbuf, Value* incount,
                                       Value* outbuf) {
    codegenPrimitive(inbuf, incount, outbuf, this->size, this->type);
}

Datatype *PrimitiveDatatype::compress() {
    return new PrimitiveDatatype(this->type);
}

void PrimitiveDatatype::globalCodegen(llvm::Module *module) {

};

int PrimitiveDatatype::getExtent() {
    return this->extent;
}

int PrimitiveDatatype::getSize() {
    return this->size;
}

string PrimitiveDatatype::toString(bool summary) {
    string res = "";
    switch (this->type) {
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


/* Class ContiguousDatatype */
ContiguousDatatype::ContiguousDatatype(int count, Datatype* basetype) {
    this->basetype = basetype->clone();
    this->count = count;
}

ContiguousDatatype::~ContiguousDatatype(void) {
    delete this->basetype;
}

ContiguousDatatype* ContiguousDatatype::clone() {
    return new ContiguousDatatype(this->count, this->basetype);
}

void ContiguousDatatype::packCodegen(Value* inbuf, Value* incount, Value* outbuf) {
    codegenContiguous(inbuf, incount, outbuf, this->basetype, 
                      this->basetype->getExtent(), this->basetype->getSize(), 
                      this->count, true);
}

void ContiguousDatatype::unpackCodegen(Value* inbuf, Value* incount, Value* outbuf) {
    codegenContiguous(inbuf, incount, outbuf, this->basetype, 
                      this->basetype->getSize(), this->basetype->getExtent(), 
                      this->count, false);
}

Datatype *ContiguousDatatype::compress() {
    Datatype *cbasetype = this->basetype->compress();
    Datatype *datatype;

    // Compress contiguous basetypes into count.  Only applies if
    // basetype's extent is the same as its size.
    ContiguousDatatype *ctg = dynamic_cast<ContiguousDatatype*>(cbasetype);
    if (ctg != NULL && ctg->getExtent() == ctg->getSize()) {
        datatype = new ContiguousDatatype(this->count * ctg->getCount(), ctg->getBasetype());
        delete cbasetype;
    }
    else {
        datatype = new ContiguousDatatype(this->count, cbasetype);
    }

    assert(datatype->getSize() == this->getSize());
    assert(datatype->getExtent() == this->getExtent());

    return datatype;
}

void ContiguousDatatype::globalCodegen(llvm::Module *mod) {
    basetype->globalCodegen(mod);
};

int ContiguousDatatype::getExtent() {
    return this->count * this->basetype->getExtent();
}

int ContiguousDatatype::getSize() {
    return this->count * this->basetype->getSize();
}

int ContiguousDatatype::getCount() {
    return this->count;
}

Datatype *ContiguousDatatype::getBasetype() {
    return this->basetype;
}

string ContiguousDatatype::toString(bool summary) {
    stringstream res;
    res << "ctg(" << this->count << ")[" << basetype->toString(summary) << "]";
    return res.str();
}


/* Class VectorDatatype */
VectorDatatype::VectorDatatype(int count, int blocklen, int stride, Datatype* basetype) {
    this->count = count;
    this->blocklen = blocklen;
    this->stride = stride;
    this->basetype = basetype->clone();
}

VectorDatatype::~VectorDatatype(void) {
    delete this->basetype;
}

VectorDatatype* VectorDatatype::clone() {
    return new VectorDatatype(this->count, this->blocklen,
                              this->stride, this->basetype);
}

void VectorDatatype::packCodegen(Value* inbuf, Value* incount, Value* outbuf) {
    codegenVector(inbuf, incount, outbuf, this->basetype, this->count,
                  this->blocklen, this->basetype->getExtent() * this->stride,
                  this->basetype->getSize() * this->blocklen, true);
}

void VectorDatatype::unpackCodegen(Value* inbuf, Value* incount, Value* outbuf) {
    codegenVector(inbuf, incount, outbuf, this->basetype, this->count,
                  this->blocklen, this->basetype->getSize() * this->blocklen,
                  this->basetype->getExtent() * this->stride, false);
}

Datatype *VectorDatatype::compress() {
    Datatype *cbasetype = this->basetype->compress();
    Datatype *datatype;

    // Compress contiguous basetypes into blocklen.  Only applies if
    // basetype's extent is the same as its size.
    ContiguousDatatype *ctg = dynamic_cast<ContiguousDatatype*>(cbasetype);
    if (ctg != NULL && ctg->getExtent() == ctg->getSize()) {
        datatype = new VectorDatatype(this->count,
                                      this->blocklen * ctg->getCount(),
                                      this->stride * ctg->getCount(),
                                      ctg->getBasetype());
        delete cbasetype;
    }
    else {
        datatype = new VectorDatatype(this->count,
                                      this->blocklen,
                                      this->stride,
                                      cbasetype);
    }

    // Try to promote the vector to a contiguous type
    VectorDatatype *vecdt = (VectorDatatype*)datatype;
    if (vecdt->getStride() == vecdt->getBlocklen()) {
        int contigCount = vecdt->getCount() * vecdt->getBlocklen();
        datatype = new ContiguousDatatype(contigCount, vecdt->getBasetype());
        delete vecdt;
    }

    assert(datatype->getSize() == this->getSize());
    assert(datatype->getExtent() == this->getExtent());

    return datatype;
}

void VectorDatatype::globalCodegen(llvm::Module *mod) {
    basetype->globalCodegen(mod);
};

int VectorDatatype::getExtent() {
    return (this->count - 1) * this->basetype->getExtent() *
        this->stride + this->blocklen * this->basetype->getExtent();
}

int VectorDatatype::getSize() {
    return this->count * this->blocklen*this->basetype->getSize();
}

int VectorDatatype::getCount() {
    return count;
}

int VectorDatatype::getBlocklen() {
    return blocklen;
}

int VectorDatatype::getStride() {
    return stride;
}

Datatype *VectorDatatype::getBasetype() {
    return this->basetype;
}

string VectorDatatype::toString(bool summary) {
    stringstream res;
    string sep = (summary) ? "," : " ";
    res << "vec(" << this->count << sep << this->blocklen << sep
        << this->stride << ")[" << basetype->toString(summary) << "]";
    return res.str();
}


/* Class HVectorDatatype */
HVectorDatatype::HVectorDatatype(int count, int blocklen, int stride, Datatype* basetype) {
    this->count = count;
    this->blocklen = blocklen;
    this->stride = stride;
    this->basetype = basetype->clone();
} 

HVectorDatatype* HVectorDatatype::clone() {
    return new HVectorDatatype(this->count, this->blocklen, this->stride, this->basetype);
}

HVectorDatatype::~HVectorDatatype(void) {
    delete this->basetype;
}

void HVectorDatatype::packCodegen(Value* inbuf, Value* incount, Value* outbuf) {
    codegenVector(inbuf, incount, outbuf, this->basetype, this->count, this->blocklen,
            this->stride, this->basetype->getSize() * this->blocklen, true);
}

void HVectorDatatype::unpackCodegen(Value* inbuf, Value* incount, Value* outbuf) {
    codegenVector(inbuf, incount, outbuf, this->basetype, this->count, this->blocklen,
            this->basetype->getSize() * this->blocklen, this->stride, false);
}

Datatype *HVectorDatatype::compress() {
    Datatype *cbasetype = this->basetype->compress();
    Datatype *datatype;

    // Compress contiguous basetypes into blocklen.  Only applies if
    // basetype's extent is the same as its size.
    ContiguousDatatype *ctg = dynamic_cast<ContiguousDatatype*>(cbasetype);
    if (ctg != NULL && ctg->getExtent() == ctg->getSize()) {
        datatype = new HVectorDatatype(this->count,
                                       this->blocklen * ctg->getCount(),
                                       this->stride,
                                       ctg->getBasetype());
        delete cbasetype;
    }
    else {
        datatype = new HVectorDatatype(this->count,
                                       this->blocklen,
                                       this->stride,
                                       cbasetype);
    }

    // Try to promote the vector to a contiguous type
    HVectorDatatype *vecdt = (HVectorDatatype*)datatype;
    if (vecdt->getStride() == (vecdt->getBlocklen() * vecdt->getBasetype()->getExtent())) {
        int contigCount = vecdt->getCount() * vecdt->getBlocklen();
        datatype = new ContiguousDatatype(contigCount, vecdt->getBasetype());
        delete vecdt;
    }

    assert(datatype->getSize() == this->getSize());
    assert(datatype->getExtent() == this->getExtent());

    return datatype;
}

void HVectorDatatype::globalCodegen(llvm::Module *mod) {
    basetype->globalCodegen(mod);
};

int HVectorDatatype::getExtent() {
    if (this->stride > 0) return (this->count-1)*this->stride + this->blocklen*this->basetype->getExtent();
    else return (-((this->count-1)*this->stride + this->blocklen*this->basetype->getExtent()) + this->count*this->blocklen*this->basetype->getExtent());
}

int HVectorDatatype::getSize() {
    return this->count * this->blocklen*this->basetype->getSize();
}

int HVectorDatatype::getCount() {
    return count;
}

int HVectorDatatype::getBlocklen() {
    return blocklen;
}

int HVectorDatatype::getStride() {
    return stride;
}

Datatype *HVectorDatatype::getBasetype() {
    return this->basetype;
}

string HVectorDatatype::toString(bool summary) {
    stringstream res;
    string sep = (summary) ? "," : " ";
    res << "hvec(" << this->count << sep << this->blocklen << sep << this->stride << ")[" << basetype->toString(summary) << "]";
    return res.str();
}


/* Class IndexedBlockDatatype */
IndexedBlockDatatype::IndexedBlockDatatype(int count, int blocklen, int* displ, Datatype* basetype) : Datatype() {
    this->count = count;
    this->basetype = basetype->clone();
    this->blocklen = blocklen;
    for (int i=0; i<count; i++) this->displs.push_back(displ[i]);
    indices_arr = NULL;
}

IndexedBlockDatatype::~IndexedBlockDatatype(void) {
    delete this->basetype;
}

IndexedBlockDatatype* IndexedBlockDatatype::clone() {
    IndexedBlockDatatype* t_new = new IndexedBlockDatatype(this->count, this->blocklen, &(this->displs[0]), this->basetype);
    return t_new;
}

void IndexedBlockDatatype::packCodegen(Value* inbuf, Value* incount, Value* outbuf) {
    codegenIndexedBlock(outbuf, inbuf, incount, this->getExtent(), this->getSize(), this->count,
                        this->blocklen, this->basetype, this->displs, indices_arr, true);
}

void IndexedBlockDatatype::unpackCodegen(Value* inbuf, Value* incount, Value* outbuf) {
    codegenIndexedBlock(inbuf, outbuf, incount, this->getExtent(), this->getSize(), this->count,
                        this->blocklen, this->basetype, this->displs, indices_arr, false);
}

Datatype *IndexedBlockDatatype::compress() {
    Datatype *cbasetype = this->basetype->compress();
    return new IndexedBlockDatatype(this->count, this->blocklen,
                                    &(this->displs[0]), cbasetype);
}

void IndexedBlockDatatype::globalCodegen(llvm::Module *mod) {
    if(count > IDXB_LOOP_TRESHOLD) {
        ArrayType* indices_types = ArrayType::get(LLVM_INT32, count);
        this->indices_arr = new GlobalVariable(*mod, indices_types, true,
                                               GlobalValue::InternalLinkage,
                                               0, "displacements");
        indices_arr->setAlignment(4);

        std::vector<Constant*> indices_vals(count);
        for (int i=0; i<count; i++) {
            indices_vals[i] = constNode(displs[i] * basetype->getExtent());
        }
        Constant* indices_initializer = ConstantArray::get(indices_types, indices_vals);
        indices_arr->setInitializer(indices_initializer);
    }

    basetype->globalCodegen(mod);
}

int IndexedBlockDatatype::getExtent() {
    if (this->count == 0) return 0;

    int bext = this->basetype->getExtent();

    int ub = this->displs[0] * bext + this->blocklen * bext;
    int lb = this->displs[0] * bext;

    for (int i=0; i<this->count; i++) {
        int tmp_ub = this->displs[i] * bext + this->blocklen * bext;
        int tmp_lb = this->displs[i] * bext;
        if (tmp_ub > ub) ub = tmp_ub;
        if (tmp_lb < lb) lb = tmp_lb;
    }

    return ub - lb;
}

int IndexedBlockDatatype::getSize() {
    int sum = 0;
    int bsize = this->basetype->getSize();
    for (int i=0; i<this->count; i++) {
        sum += bsize * this->blocklen;
    }

    return sum;
}

string IndexedBlockDatatype::toString(bool summary) {
    stringstream res;
    string sep = (summary) ? "," : " ";
    if (summary) {
        res << "idxb(" << this->count << "," << this->blocklen;
    }
    else {
        res << "idxb("  << this->blocklen << ":";
        for (unsigned int i=0; i<displs.size(); i++) {
            res << displs[i];
            if (i <displs.size() - 1) {
                res << sep;
            }
        }
    }
    res << ")[" << basetype->toString(summary) << "]";
    return res.str();
}


/* Class HIndexedDatatype */
HIndexedDatatype::HIndexedDatatype(int count, int* blocklen, long* displ, Datatype* basetype) : Datatype() {
    this->count = count;
    this->basetype = basetype->clone();
    for (int i=0; i<count; i++) this->blocklens.push_back(blocklen[i]);
    for (int i=0; i<count; i++) this->displs.push_back(displ[i]);
}

HIndexedDatatype::~HIndexedDatatype(void) {
    delete this->basetype;
}

HIndexedDatatype* HIndexedDatatype::clone() {
    HIndexedDatatype* t_new = new HIndexedDatatype(this->count, &(this->blocklens[0]), &(this->displs[0]), this->basetype);
    return t_new;
}

void HIndexedDatatype::packCodegen(Value* inbuf, Value* incount, Value* outbuf) {
    codegenHindexed(outbuf, inbuf, incount, this->getExtent(), this->count,
                    this->basetype, this->blocklens, this->displs, true);
}

void HIndexedDatatype::unpackCodegen(Value* inbuf, Value* incount, Value* outbuf) {
    codegenHindexed(inbuf, outbuf, incount, this->getExtent(), this->count,
                    this->basetype, this->blocklens, this->displs, false);
}

Datatype *HIndexedDatatype::compress() {
    Datatype *cbasetype = this->basetype->compress();
    return new HIndexedDatatype(this->count, &(this->blocklens[0]),
                                &(this->displs[0]), cbasetype);
}

void HIndexedDatatype::globalCodegen(llvm::Module *mod) {
    basetype->globalCodegen(mod);
};

int HIndexedDatatype::getExtent() {
    if (this->count == 0) return 0;

    int bext = this->basetype->getExtent();
    int ub = this->displs[0] + bext * this->blocklens[0];
    int lb = this->displs[0];

    for (int i=0; i<this->count; i++) {
        int tmp_ub = this->displs[i] + bext * this->blocklens[i];
        int tmp_lb = this->displs[i];
        if (tmp_ub > ub) ub = tmp_ub;
        if (tmp_lb < lb) lb = tmp_lb;
    }

    return ub - lb;
}

int HIndexedDatatype::getSize() {
    int sum = 0;
    int bsize = this->basetype->getSize();
    for (int i=0; i<this->count; i++) {
        sum += bsize * this->blocklens[i];
    }

    return sum;
}

string HIndexedDatatype::toString(bool summary) {
    stringstream res;
    string sep = (summary) ? ";" : " ";
    res << "hidx(";
    for (unsigned int i=0; i<displs.size(); i++) {
        res << displs[i] << "," << blocklens[i];
        if (i <displs.size() - 1) {
            res << sep;
        }
    }
    res << ")[" << basetype->toString(summary) << "]";
    return res.str();
}


/* Class StructDatatype */
StructDatatype::StructDatatype(int count, int* blocklen, long*  displ,
                               Datatype** types) : Datatype() {
    this->count = count;
    for (int i=0; i<count; i++) blocklens.push_back(blocklen[i]);
    for (int i=0; i<count; i++) displs.push_back(displ[i]);
    for (int i=0; i<count; i++) basetypes.push_back(types[i]->clone());
}

StructDatatype::~StructDatatype(void) {
    for (int i=0; i<this->count; i++) delete basetypes[i];
}

StructDatatype* StructDatatype::clone() {
    StructDatatype* t_new = new StructDatatype(this->count, &(this->blocklens[0]), &(this->displs[0]), &(this->basetypes[0]));
    return t_new;
}

void StructDatatype::packCodegen(Value* inbuf, Value* incount, Value* outbuf) {
    codegenStruct(outbuf, inbuf, incount, this->getExtent(), this->count,
                  this->blocklens, this->displs, this->basetypes, true);
}

void StructDatatype::unpackCodegen(Value* inbuf, Value* incount, Value* outbuf) {
    codegenStruct(inbuf, outbuf, incount, this->getExtent(), this->count,
                  this->blocklens, this->displs, this->basetypes, false);
}

Datatype *StructDatatype::compress() {
    std::vector<Datatype*> cbasetypes(this->basetypes.size());
    for (unsigned int i=0 ; i<this->basetypes.size(); i++) {
        cbasetypes[i] = this->basetypes[i]->compress();
    }
    return new StructDatatype(this->count, &(this->blocklens[0]),
                                  &(this->displs[0]), &(this->basetypes[0]));
}

void StructDatatype::globalCodegen(llvm::Module *mod) {
    for (unsigned int i=0 ; i<this->basetypes.size(); i++) {
        basetypes[i]->globalCodegen(mod);
    }
};

int StructDatatype::getExtent() {
    if (this->count == 0) return 0;

    int lb = this->displs[0];
    int ub = this->displs[0] + this->basetypes[0]->getExtent() * this->blocklens[0];
    for (int i=0; i<this->count; i++) {
        int tmp_ub = this->displs[i] + this->basetypes[i]->getExtent() * this->blocklens[i];
        int tmp_lb = this->displs[i];
        if (tmp_ub > ub) ub = tmp_ub;
        if (tmp_lb < lb) lb = tmp_lb;
    }

    return ub-lb;
}

int StructDatatype::getSize() {
    int sum = 0;
    for (int i=0; i<this->count; i++) {
        sum += this->basetypes[i]->getSize() * this->blocklens[i];
    }

    return sum;
}

string StructDatatype::toString(bool summary) {
    stringstream res;
    string sep = (summary) ? ";" : " ";
    res << "struct(" << this->count << ")";
    for (unsigned int i=0; i<displs.size(); i++) {
        res << "[" << this->displs[i] << "," << this->blocklens[i] << "]";
        res << "{" << this->basetypes[i]->toString(summary) << "}";
        if (i <displs.size() - 1) {
            res << sep;
        }
    }
    return res.str();
}

void DDT_Commit(Datatype* ddt) {
#if DDT_OUTPUT
    ddt->print();
#endif
#if !LAZY
    ddt->compile(Datatype::PACK_UNPACK);
#endif
}

// this calls the pack/unpack function
void DDT_Pack(void* inbuf, void* outbuf, Datatype* ddt, int count) {
#if LAZY
    if (ddt->pack == NULL) ddt->compile(Datatype::PACK);
#endif
    ddt->pack(inbuf, count, outbuf);
}

void DDT_Lazy_Unpack_Commit(Datatype* ddt) {
#if LAZY
    if (ddt->unpack == NULL) ddt->compile(Datatype::UNPACK);
#endif
}

void DDT_Unpack(void* inbuf, void* outbuf, Datatype* ddt, int count) {
#if LAZY
    DDT_Lazy_Unpack_Commit(ddt);
#endif
    ddt->unpack(inbuf, count, outbuf);
}

void DDT_Free(Datatype* ddt) {
    delete ddt;
}

// init the JIT compiler
void DDT_Init() {
    InitializeNativeTarget();
    LLVMContext &Context = getGlobalContext();
    module = new Module("FARC-JIT", Context);

    // Create the JIT.  This takes ownership of the module.
    std::string ErrStr;
    EngineBuilder engine_builder(module);
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
    FunctionPassManager* OurFPM = new FunctionPassManager(module);

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

}

} // namespace farc

/* Define the functions of the pack header file */
