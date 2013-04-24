// Copyright 2013 Timo Schneider and Fredrik Berg Kjolstad
//
// This file is part of the libpack packing library.
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT in the top level directory for details.

#ifndef DDT_JIT_H
#define DDT_JIT_H

#include <cstdlib>
#include <vector>
#include <string>

/* Forward declare llvm values */
namespace llvm {
class Value;
class Function;
class Module;
class GlobalVariable;
}

namespace farc {

enum DatatypeName {PRIMITIVE, CONTIGUOUS, VECTOR, HVECTOR, INDEXEDBLOCK, HINDEXED, STRUCT, RESIZED};

/* Base class for all datatypes */
class Datatype {
public:
    enum CompilationType { PACK, UNPACK, PACK_UNPACK };

    Datatype() { this->pack = NULL; this->unpack = NULL; }
    virtual ~Datatype();
    virtual Datatype* clone() = 0;

    virtual int getExtent() = 0;
    virtual int getTrueExtent() = 0;
    virtual int getSize() = 0;
    virtual int getLowerBound() = 0;
    virtual int getTrueLowerBound() = 0;
    virtual int getUpperBound() = 0;
    virtual int getTrueUpperBound() = 0;
    virtual DatatypeName getDatatypeName() = 0;
    virtual std::vector<Datatype*> getSubtypes() = 0;


	virtual std::string toString(bool summary = false) = 0;
    virtual void print(bool summary = false);

    virtual void compile(CompilationType type);
    void (*pack)(void*, int, void*);
    void (*unpack)(void*, int, void*);

    virtual Datatype *compress() = 0;
    virtual void packCodegen(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf) = 0;
    virtual void unpackCodegen(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf) = 0;
    virtual void globalCodegen(llvm::Module *mod) = 0;

    // TODO: remove this function
    virtual void cleanup() {}

    llvm::Function* fpack;
    llvm::Function* funpack;
};

/* Class for primitive types, such as MPI_INT, MPI_BYTE, etc */
class PrimitiveDatatype : public Datatype {
public:
    enum PrimitiveType { BYTE, CHAR, DOUBLE, FLOAT, INT };

    PrimitiveDatatype(PrimitiveType type);
    virtual ~PrimitiveDatatype(void) {};
    PrimitiveDatatype* clone();

    std::vector<Datatype*> getSubtypes();
    DatatypeName getDatatypeName();
    int getExtent();
    int getTrueExtent();
    int getSize();
    int getLowerBound();
    int getTrueLowerBound();
    int getUpperBound();
    int getTrueUpperBound();

    std::string toString(bool summary = false);

    Datatype *compress();
    void packCodegen(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void unpackCodegen(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void globalCodegen(llvm::Module *mod);

private:
    int size;
    int lower_bound;
    int upper_bound;
    int true_lower_bound;
    int true_upper_bound;

    PrimitiveDatatype::PrimitiveType type;
};

/* Class for contiguous types */
class ContiguousDatatype : public Datatype {
public:
    ContiguousDatatype(int count, Datatype* basetype);
    virtual ~ContiguousDatatype(void);
    ContiguousDatatype* clone();

    std::vector<Datatype*> getSubtypes();
    DatatypeName getDatatypeName();
    int getExtent();
    int getTrueExtent();
    int getSize();
    int getLowerBound();
    int getTrueLowerBound();
    int getUpperBound();
    int getTrueUpperBound();

    int getCount();
    Datatype *getBasetype();
    std::string toString(bool summary = false);

    Datatype *compress();
    void packCodegen(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void unpackCodegen(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void globalCodegen(llvm::Module *mod);

private:
    int size;
    int lower_bound;
    int upper_bound;
    int true_lower_bound;
    int true_upper_bound;

    int count;
    Datatype* basetype;
};

/* Class for vector types */
class VectorDatatype : public Datatype {
public:
    VectorDatatype(int count, int blocklen, int stride, Datatype* basetype);
    virtual ~VectorDatatype(void);
    VectorDatatype* clone();

    std::vector<Datatype*> getSubtypes();
    DatatypeName getDatatypeName();
    int getExtent();
    int getTrueExtent();
    int getSize();
    int getLowerBound();
    int getTrueLowerBound();
    int getUpperBound();
    int getTrueUpperBound();

    int getCount();
    int getBlocklen();
    int getStride();
    Datatype *getBasetype();
    std::string toString(bool summary = false);

    Datatype *compress();
    void packCodegen(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void unpackCodegen(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void globalCodegen(llvm::Module *mod);

private:
    int size;
    int lower_bound;
    int upper_bound;
    int true_lower_bound;
    int true_upper_bound;

    int count;
    int blocklen;
    int stride;
    Datatype* basetype;
};

/* Class for hvector types */
class HVectorDatatype : public Datatype {
public:
    HVectorDatatype(int count, int blocklen, int stride, Datatype* basetype);
    virtual ~HVectorDatatype(void);
    HVectorDatatype* clone();

    std::vector<Datatype*> getSubtypes();
    DatatypeName getDatatypeName();
    int getExtent();
    int getTrueExtent();
    int getSize();
    int getLowerBound();
    int getTrueLowerBound();
    int getUpperBound();
    int getTrueUpperBound();

    int getCount();
    int getBlocklen();
    int getStride();
    Datatype *getBasetype();
    std::string toString(bool summary = false);

    Datatype *compress();
    void packCodegen(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void unpackCodegen(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void globalCodegen(llvm::Module *mod);

private:
    int size;
    int lower_bound;
    int upper_bound;
    int true_lower_bound;
    int true_upper_bound;

    int count;
    int blocklen;
    int stride;
    Datatype* basetype;
};

/* Class for indexed block types */
class IndexedBlockDatatype : public Datatype {
public:
    IndexedBlockDatatype(int count, int blocklen, int* displ, Datatype* basetype);
    virtual ~IndexedBlockDatatype(void);
    IndexedBlockDatatype* clone();

    std::vector<Datatype*> getSubtypes();
    DatatypeName getDatatypeName();
    int getExtent();
    int getTrueExtent();
    int getSize();
    int getLowerBound();
    int getTrueLowerBound();
    int getUpperBound();
    int getTrueUpperBound();

    std::string toString(bool summary = false);

    Datatype *compress();
    void packCodegen(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void unpackCodegen(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void globalCodegen(llvm::Module *mod);

private:
    int size;
    int lower_bound;
    int upper_bound;
    int true_lower_bound;
    int true_upper_bound;

    int count;
    int blocklen;
    std::vector<int> displs;
    Datatype* basetype;

    llvm::GlobalVariable* indices_arr;

    // TODO Remove this funciton
    void cleanup();
};

/* Class for hindexed types */
class HIndexedDatatype : public Datatype {
public:
    HIndexedDatatype(int count, int* blocklen, long* displ, Datatype* basetype);
    virtual ~HIndexedDatatype(void);
    HIndexedDatatype* clone();

    std::vector<Datatype*> getSubtypes();
    DatatypeName getDatatypeName();
    int getExtent();
    int getTrueExtent();
    int getSize();
    int getLowerBound();
    int getTrueLowerBound();
    int getUpperBound();
    int getTrueUpperBound();

    std::string toString(bool summary = false);

    Datatype *compress();
    void packCodegen(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void unpackCodegen(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void globalCodegen(llvm::Module *mod);

private:
    int size;
    int lower_bound;
    int upper_bound;
    int true_lower_bound;
    int true_upper_bound;

    int count;
    std::vector<int> blocklens;
    std::vector<long> displs;
    Datatype* basetype;
};

/* Class for struct types */
class StructDatatype : public Datatype {
public:
    StructDatatype(int count, int* blocklen, long*  displ, Datatype** basetypes);
    virtual ~StructDatatype(void);
    StructDatatype* clone();

    std::vector<Datatype*> getSubtypes();
    DatatypeName getDatatypeName();
    int getExtent();
    int getTrueExtent();
    int getSize();
    int getLowerBound();
    int getTrueLowerBound();
    int getUpperBound();
    int getTrueUpperBound();

    std::string toString(bool summary = false);

    Datatype *compress();
    void packCodegen(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void unpackCodegen(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void globalCodegen(llvm::Module *mod);

private:
    int size;
    int lower_bound;
    int upper_bound;
    int true_lower_bound;
    int true_upper_bound;

    int count;
    std::vector<int> blocklens;
    std::vector<long> displs;
    std::vector<Datatype*> basetypes;
};

/* Class for resized types */
class ResizedDatatype : public Datatype {
public:
    ResizedDatatype(Datatype* basetype, int lb, int extent);
    virtual ~ResizedDatatype(void);
    ResizedDatatype* clone();

    std::vector<Datatype*> getSubtypes();
    DatatypeName getDatatypeName();
    int getExtent();
    int getTrueExtent();
    int getSize();
    int getLowerBound();
    int getTrueLowerBound();
    int getUpperBound();
    int getTrueUpperBound();

    std::string toString(bool summary = false);

    Datatype *compress();
    void packCodegen(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void unpackCodegen(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void globalCodegen(llvm::Module *mod);

private:
    int size;
    int lower_bound;
    int upper_bound;
    int true_lower_bound;
    int true_upper_bound;

    Datatype* basetype;
};




/* FARC Library Functions */
void DDT_Init();
void DDT_Finalize();

void DDT_Commit(Datatype* ddt);
void DDT_Lazy_Unpack_Commit(Datatype* ddt);  // This function should be removed
void DDT_Free(Datatype* ddt);

void DDT_Pack(void* inbuf, void* outbuf, Datatype* ddt, int count);
void DDT_Unpack(void* inbuf, void* outbuf, Datatype* ddt, int count);

void DDT_Pack_partial(void* inbuf, void* outbuf, Datatype* ddt, int count, int segnum);
void DDT_Unpack_partial(void* inbuf, void* outbuf, Datatype* ddt, int count, int segnum);

} // namespace farc

#endif

