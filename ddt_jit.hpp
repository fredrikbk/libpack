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
    virtual void packCodegen(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf) = 0;
    virtual void unpackCodegen(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf) = 0;
    void (*pack)(void*, int, void*);
    void (*unpack)(void*, int, void*);

    llvm::Function* FPack;
    llvm::Function* FUnpack;
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

private:
    PrimitiveDatatype::PrimitiveType Type;
    int Extent;
    int Size;
};

/* Class for contiguous types */
class ContiguousDatatype : public Datatype {
public:
    ContiguousDatatype(Datatype* type, int count);
    virtual ~ContiguousDatatype(void);
    ContiguousDatatype* clone();

    int getExtent();
    int getSize();
    std::string toString();

    void packCodegen(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void unpackCodegen(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);

private:
    Datatype* Basetype;
    int Count;
};

/* Class for vector types */
class VectorDatatype : public Datatype {
public:
    VectorDatatype(Datatype* type, int count, int blocklen, int stride); 
    virtual ~VectorDatatype(void); 
    VectorDatatype* clone();

    int getExtent();
    int getSize();
    std::string toString();

    void packCodegen(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void unpackCodegen(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);

private:
    Datatype* Basetype;
    int Count;
    int Blocklen;
    int Stride;
};

/* Class for hvector types */
class HVectorDatatype : public Datatype {
public:
    HVectorDatatype(Datatype* type, int count, int blocklen, int stride);
    virtual ~HVectorDatatype(void);
    HVectorDatatype* clone();

    int getExtent();
    int getSize();
    std::string toString();

    void packCodegen(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void unpackCodegen(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);

private:
    Datatype* Basetype;
    int Count;
    int Blocklen;
    int Stride;
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

private:
    int Count;
    Datatype* Basetype;
    int Blocklen;
    std::vector<int> Displ;
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

private:
    int Count;
    Datatype* Basetype;
    std::vector<int> Blocklen;
    std::vector<long> Displ;
};

/* Class for struct types */
class StructDatatype : public Datatype {
public:
    StructDatatype(int count, int* blocklen, long*  displ, Datatype** types);
    virtual ~StructDatatype(void);
    StructDatatype* clone();

    int getExtent();
    int getSize();
    std::string toString();

    void packCodegen(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void unpackCodegen(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);

private:
    int Count;
    std::vector<Datatype*> Types;
    std::vector<int> Blocklen;
    std::vector<long> Displ;
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

