#ifndef DDT_JIT_H
#define DDT_JIT_H

#include <vector>
#include <mpi.h>

#if ((__clang_major__ == 3) && (__clang_minor__ == 2))
#define LLVM32 1
#warn foo
#else
#define LLVM32 0
#endif

#if LLVM32
#include <llvm/Value.h>
#else
#include <llvm/IR/Value.h>
#endif

#define LAZY 0 
#define TIME 0 

/* Base class for all datatypes */
class FARC_Datatype {

    public:
    FARC_Datatype() { this->packer = NULL; this->unpacker = NULL; }
    virtual ~FARC_Datatype() {}
    virtual llvm::Value *Codegen_Pack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf) = 0;
    virtual llvm::Value *Codegen_Unpack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf) = 0;
    virtual int getExtend() = 0;
    virtual int getSize() = 0;

    void (*packer)(void*, int, void*);
    void (*unpacker)(void*, int, void*);

};

/* Class for primitive types, such as MPI_INT, MPI_BYTE, etc */
class FARC_PrimitiveDatatype : public FARC_Datatype {

    MPI_Datatype Type; // this MUST be a primitive type
    int Extend;
    int Size;

    public:
    FARC_PrimitiveDatatype(MPI_Datatype type);
    llvm::Value* Codegen_Pack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    llvm::Value* Codegen_Unpack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    int getExtend();
    int getSize();

};

/* Class for contiguous types */
class FARC_ContiguousDatatype : public FARC_Datatype {

    FARC_Datatype* Basetype;
    int Count;
    void Codegen(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf, int elemstride_in, int elemstride_out, bool pack);

    public:
    FARC_ContiguousDatatype(FARC_Datatype* type, int count) : FARC_Datatype(), Basetype(type), Count(count) {}
    int getExtend();
    int getSize();
    llvm::Value *Codegen_Pack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    llvm::Value *Codegen_Unpack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);

};

/* Class for vector types */
class FARC_VectorDatatype : public FARC_Datatype {

    FARC_Datatype* Basetype;
    int Count;
    int Blocklen;
    int Stride;

    public:
    FARC_VectorDatatype(FARC_Datatype* type, int count, int blocklen, int stride) : FARC_Datatype(), Basetype(type), Count(count), Blocklen(blocklen), Stride(stride) {}
    llvm::Value *Codegen_Pack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    llvm::Value *Codegen_Unpack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    int getExtend();
    int getSize();

};

/* Class for hvector types */
class FARC_HVectorDatatype : public FARC_Datatype {

    FARC_Datatype* Basetype;
    int Count;
    int Blocklen;
    int Stride;

    public:
    FARC_HVectorDatatype(FARC_Datatype* type, int count, int blocklen, int stride) :  FARC_Datatype(), Basetype(type), Count(count), Blocklen(blocklen), Stride(stride) {}
    llvm::Value *Codegen_Pack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    llvm::Value *Codegen_Unpack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    int getExtend();
    int getSize();

};

/* Class for hindexed types */
class FARC_HIndexedDatatype : public FARC_Datatype {

    int Count;
    FARC_Datatype* Basetype;
    std::vector<int> Blocklen;
    std::vector<int> Displ;
    void Codegen(llvm::Value *contig_buf, llvm::Value *hindexed_buf, llvm::Value* incount, bool gather);

    public:
    FARC_HIndexedDatatype(int count, int* blocklen, long* displ, FARC_Datatype* basetype);
    llvm::Value *Codegen_Pack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    llvm::Value *Codegen_Unpack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    int getExtend();
    int getSize();

};

/* Class for struct types */
class FARC_StructDatatype : public FARC_Datatype {

    int Count;
    std::vector<FARC_Datatype*> Types;
    std::vector<int> Blocklen;
    std::vector<int> Displ;

    public:
    FARC_StructDatatype(int count, int* blocklen, long*  displ, FARC_Datatype** types);
    llvm::Value *Codegen_Pack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    llvm::Value *Codegen_Unpack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    int getExtend();
    int getSize();

};


/* FARC Library Functions */
void FARC_DDT_Init();
void FARC_DDT_Finalize();

FARC_Datatype* FARC_DDT_Commit(FARC_Datatype* ddt);
void FARC_DDT_Lazy_Unpack_Commit(FARC_Datatype* ddt);  // This function should be removed
void FARC_DDT_Free(FARC_Datatype* ddt);

void FARC_DDT_Pack(void* inbuf, void* outbuf, FARC_Datatype* ddt, int count);
void FARC_DDT_Unpack(void* inbuf, void* outbuf, FARC_Datatype* ddt, int count);

#endif

