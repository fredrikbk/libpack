#include "ddt_jit.hpp"
#include "pack.h"

#include "codegen.hpp"
#include "codegen_common.hpp"

#include <map>
#include <cstdio>
#include <sstream>

#include "llvm/IR/Module.h"
#include "llvm/IR/DerivedTypes.h"
#include "llvm/IR/Intrinsics.h"

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

using namespace std;
using namespace llvm;

namespace farc {

static Module *TheModule;
static std::map<std::string, Value*> NamedValues;
static ExecutionEngine *TheExecutionEngine;

std::vector<std::string> Args;
FunctionType *FT;


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

void PrimitiveDatatype::Codegen_Pack(Value* inbuf, Value* incount,
                                     Value* outbuf) {
    codegenPrimitive(inbuf, incount, outbuf, this->Size, this->Type);
}

void PrimitiveDatatype::Codegen_Unpack(Value* inbuf, Value* incount,
                                       Value* outbuf) {
    codegenPrimitive(inbuf, incount, outbuf, this->Size, this->Type);
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

void ContiguousDatatype::Codegen_Pack(Value* inbuf, Value* incount, Value* outbuf) {
    codegenContiguous(inbuf, incount, outbuf, this->Basetype, 
                      this->Basetype->getExtent(), this->Basetype->getSize(), 
                      this->Count, true);
}

void ContiguousDatatype::Codegen_Unpack(Value* inbuf, Value* incount, Value* outbuf) {
    codegenContiguous(inbuf, incount, outbuf, this->Basetype, 
                      this->Basetype->getSize(), this->Basetype->getExtent(), 
                      this->Count, false);
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

    VectorDatatype* t_new = new VectorDatatype(this->Basetype, this->Count, 
                                               this->Blocklen, this->Stride);

    return t_new;

}

void VectorDatatype::Codegen_Pack(Value* inbuf, Value* incount, Value* outbuf) {
    codegenVector(inbuf, incount, outbuf, this->Basetype, this->Count,
                  this->Blocklen, this->Basetype->getExtent() * this->Stride,
                  this->Basetype->getSize() * this->Blocklen, true);
}

void VectorDatatype::Codegen_Unpack(Value* inbuf, Value* incount, Value* outbuf) {
    codegenVector(inbuf, incount, outbuf, this->Basetype, this->Count,
                  this->Blocklen, this->Basetype->getSize() * this->Blocklen,
                  this->Basetype->getExtent() * this->Stride, false);
}

int VectorDatatype::getExtent() {
    return (this->Count - 1) * this->Basetype->getExtent() *
        this->Stride + this->Blocklen * this->Basetype->getExtent();
}

int VectorDatatype::getSize() {
    return this->Count * this->Blocklen*this->Basetype->getSize();
}

string VectorDatatype::toString() {
    stringstream res;
    res << "vec(" << this->Count << " " << this->Blocklen << " "
        << this->Stride << ")[" << Basetype->toString() << "]";
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
    codegenVector(inbuf, incount, outbuf, this->Basetype, this->Count, this->Blocklen,
            this->Stride, this->Basetype->getSize() * this->Blocklen, true);
}

void HVectorDatatype::Codegen_Unpack(Value* inbuf, Value* incount, Value* outbuf) {
    codegenVector(inbuf, incount, outbuf, this->Basetype, this->Count, this->Blocklen,
            this->Basetype->getSize() * this->Blocklen, this->Stride, false);
}

int HVectorDatatype::getExtent() {
    if (this->Stride > 0) return (this->Count-1)*this->Stride + this->Blocklen*this->Basetype->getExtent();
    else return (-((this->Count-1)*this->Stride + this->Blocklen*this->Basetype->getExtent()) + this->Count*this->Blocklen*this->Basetype->getExtent());
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

void IndexedBlockDatatype::Codegen_Pack(Value* inbuf, Value* incount, Value* outbuf) {
    codegenIndexedBlock(outbuf, inbuf, incount, this->getExtent(), this->Count,
                        this->Blocklen, this->Basetype, this->Displ, true);
}

void IndexedBlockDatatype::Codegen_Unpack(Value* inbuf, Value* incount, Value* outbuf) {
    codegenIndexedBlock(inbuf, outbuf, incount, this->getExtent(), this->Count,
                        this->Blocklen, this->Basetype, this->Displ, false);
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

void HIndexedDatatype::Codegen_Pack(Value* inbuf, Value* incount, Value* outbuf) {
    codegenHindexed(outbuf, inbuf, incount, this->getExtent(), this->Count,
                    this->Basetype, this->Blocklen, this->Displ, true);
}

void HIndexedDatatype::Codegen_Unpack(Value* inbuf, Value* incount, Value* outbuf) {
    codegenHindexed(inbuf, outbuf, incount, this->getExtent(), this->Count,
                    this->Basetype, this->Blocklen, this->Displ, false);
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
StructDatatype::StructDatatype(int count, int* blocklen, long*  displ,
                               Datatype** types) : Datatype() {
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

void StructDatatype::Codegen_Pack(Value* inbuf, Value* incount, Value* outbuf) {
    codegenStruct(outbuf, inbuf, incount, this->getExtent(), this->Count,
                  this->Blocklen, this->Displ, this->Types, true);
}

void StructDatatype::Codegen_Unpack(Value* inbuf, Value* incount, Value* outbuf) {
    codegenStruct(inbuf, outbuf, incount, this->getExtent(), this->Count,
                  this->Blocklen, this->Displ, this->Types, false);
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

}

void generate_unpack_function(Datatype* ddt) {
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
}

void DDT_Commit(Datatype* ddt) {
#if DDT_OUTPUT
    ddt->print();
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

}

} // namespace farc

/* Define the functions of the pack header file */
