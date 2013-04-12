#ifndef DDT_JIT_H
#define DDT_JIT_H

#include <cstdlib>
#include <vector>
#include <string>

/* Forward declare llvm values */
namespace llvm {
    class Value;
    class Function;
}

namespace farc {

/* Base class for all datatypes */
class Datatype {
public:
    enum CompilationType { PACK, UNPACK, PACK_UNPACK };

    Datatype() { this->pack = NULL; this->unpack = NULL; }
    virtual ~Datatype();
    virtual Datatype* clone() = 0;

    virtual int getExtent() = 0;
    virtual int getSize() = 0;
	virtual std::string toString() = 0;
    virtual void print();

    virtual void compile(CompilationType type);
    void (*pack)(void*, int, void*);
    void (*unpack)(void*, int, void*);

    virtual void packCodegen(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf) = 0;
    virtual void unpackCodegen(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf) = 0;
    virtual Datatype *compress() = 0;

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

    int getExtent();
    int getSize();
    std::string toString();

    void packCodegen(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void unpackCodegen(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    Datatype *compress();

private:
    PrimitiveDatatype::PrimitiveType type;
    int extent;
    int size;
};

/* Class for contiguous types */
class ContiguousDatatype : public Datatype {
public:
    ContiguousDatatype(int count, Datatype* basetype);
    virtual ~ContiguousDatatype(void);
    ContiguousDatatype* clone();

    int getExtent();
    int getSize();
    int getCount();
    Datatype *getBasetype();
    std::string toString();

    void packCodegen(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void unpackCodegen(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    Datatype *compress();

private:
    int count;
    Datatype* basetype;
};

/* Class for vector types */
class VectorDatatype : public Datatype {
public:
    VectorDatatype(int count, int blocklen, int stride, Datatype* basetype);
    virtual ~VectorDatatype(void);
    VectorDatatype* clone();

    int getExtent();
    int getSize();
    int getCount();
    int getBlocklen();
    int getStride();
    Datatype *getBasetype();
    std::string toString();

    void packCodegen(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void unpackCodegen(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    Datatype *compress();

private:
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

    int getExtent();
    int getSize();
    int getCount();
    int getBlocklen();
    int getStride();
    Datatype *getBasetype();
    std::string toString();

    void packCodegen(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void unpackCodegen(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    Datatype *compress();

private:
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

    int getExtent();
    int getSize();
    std::string toString();

    void packCodegen(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void unpackCodegen(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    Datatype *compress();

private:
    int count;
    int blocklen;
    std::vector<int> displs;
    Datatype* basetype;
};

/* Class for hindexed types */
class HIndexedDatatype : public Datatype {
public:
    HIndexedDatatype(int count, int* blocklen, long* displ, Datatype* basetype);
    virtual ~HIndexedDatatype(void);
    HIndexedDatatype* clone();

    int getExtent();
    int getSize();
    std::string toString();

    void packCodegen(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void unpackCodegen(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    Datatype *compress();

private:
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

    int getExtent();
    int getSize();
    std::string toString();

    void packCodegen(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void unpackCodegen(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    Datatype *compress();

private:
    int count;
    std::vector<int> blocklens;
    std::vector<long> displs;
    std::vector<Datatype*> basetypes;
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

