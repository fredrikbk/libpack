#ifndef DDT_JIT_H
#define DDT_JIT_H

#include <cstdlib>
#include <vector>
#include <string>

#define LAZY 0 
#define TIME 0 

#define DDT_OUTPUT     0 
#define LLVM_OUTPUT    0 
#define LLVM_OPTIMIZE  0 

/* Forward declare llvm values */
namespace llvm {
    class Value;
    class Function;
}

/* Base class for all datatypes */
class FARC_Datatype {

    public:
    FARC_Datatype() { this->packer = NULL; this->unpacker = NULL; }
    virtual ~FARC_Datatype() {}
    virtual void Codegen_Pack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf) = 0;
    virtual void Codegen_Unpack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf) = 0;
    virtual int getExtend() = 0;
    virtual int getSize() = 0;
    virtual void print(std::string indent) = 0;

    void (*packer)(void*, int, void*);
    void (*unpacker)(void*, int, void*);

    llvm::Function* FPack;
    llvm::Function* FUnpack;
};

/* Class for primitive types, such as MPI_INT, MPI_BYTE, etc */
class FARC_PrimitiveDatatype : public FARC_Datatype {
    public:
    enum PrimitiveType { BYTE, CHAR, DOUBLE, INT };   

    private:
    FARC_PrimitiveDatatype::PrimitiveType Type;
    int Extend;
    int Size;

    public:
    FARC_PrimitiveDatatype(PrimitiveType type);
    void Codegen_Pack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void Codegen_Unpack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    int getExtend();
    int getSize();
    void print(std::string indent);

};

/* Class for contiguous types */
class FARC_ContiguousDatatype : public FARC_Datatype {

    FARC_Datatype* Basetype;
    int Count;
    void Codegen(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf, int elemstride_in, int elemstride_out, bool pack);

    public:
    FARC_ContiguousDatatype(FARC_Datatype* type, int count) : FARC_Datatype(), Basetype(type), Count(count) {}
    void Codegen_Pack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void Codegen_Unpack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    int getExtend();
    int getSize();
    void print(std::string indent);

};

/* Class for vector types */
class FARC_VectorDatatype : public FARC_Datatype {

    FARC_Datatype* Basetype;
    int Count;
    int Blocklen;
    int Stride;

    public:
    FARC_VectorDatatype(FARC_Datatype* type, int count, int blocklen, int stride) : FARC_Datatype(), Basetype(type), Count(count), Blocklen(blocklen), Stride(stride) {}
    void Codegen_Pack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void Codegen_Unpack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    int getExtend();
    int getSize();
    void print(std::string indent);

};

/* Class for hvector types */
class FARC_HVectorDatatype : public FARC_Datatype {

    FARC_Datatype* Basetype;
    int Count;
    int Blocklen;
    int Stride;

    public:
    FARC_HVectorDatatype(FARC_Datatype* type, int count, int blocklen, int stride) :  FARC_Datatype(), Basetype(type), Count(count), Blocklen(blocklen), Stride(stride) {}
    void Codegen_Pack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void Codegen_Unpack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    int getExtend();
    int getSize();
    void print(std::string indent);

};

/* Class for indexed block types */
class FARC_IndexedBlockDatatype : public FARC_Datatype {

    int Count;
    FARC_Datatype* Basetype;
    int Blocklen;
    std::vector<int> Displ;
    void Codegen(llvm::Value *contig_buf, llvm::Value *hindexed_buf, llvm::Value* incount, bool gather);

    public:
    FARC_IndexedBlockDatatype(int count, int blocklen, int* displ, FARC_Datatype* basetype);
    void Codegen_Pack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void Codegen_Unpack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    int getExtend();
    int getSize();
    void print(std::string indent);

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
    void Codegen_Pack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void Codegen_Unpack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    int getExtend();
    int getSize();
    void print(std::string indent);

};

/* Class for struct types */
class FARC_StructDatatype : public FARC_Datatype {

    int Count;
    std::vector<FARC_Datatype*> Types;
    std::vector<int> Blocklen;
    std::vector<int> Displ;

    public:
    FARC_StructDatatype(int count, int* blocklen, long*  displ, FARC_Datatype** types);
    void Codegen_Pack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void Codegen_Unpack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void Codegen(llvm::Value *compactbuf, llvm::Value *scatteredbuf, llvm::Value* incount, bool pack);
    int getExtend();
    int getSize();
    void print(std::string indent);

};


/* FARC Library Functions */
void FARC_DDT_Init();
void FARC_DDT_Finalize();

void FARC_DDT_Commit(FARC_Datatype* ddt);
void FARC_DDT_Lazy_Unpack_Commit(FARC_Datatype* ddt);  // This function should be removed
void FARC_DDT_Free(FARC_Datatype* ddt);

void FARC_DDT_Pack(void* inbuf, void* outbuf, FARC_Datatype* ddt, int count);
void FARC_DDT_Unpack(void* inbuf, void* outbuf, FARC_Datatype* ddt, int count);

void FARC_DDT_Pack_partial(void* inbuf, void* outbuf, FARC_Datatype* ddt, int count, int segnum);
void FARC_DDT_Unpack_partial(void* inbuf, void* outbuf, FARC_Datatype* ddt, int count, int segnum);


#endif

