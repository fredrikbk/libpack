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

namespace farc {

/* Base class for all datatypes */
class Datatype {

    public:
    Datatype() { this->packer = NULL; this->unpacker = NULL; }
    virtual ~Datatype() {}
    virtual void Codegen_Pack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf) = 0;
    virtual void Codegen_Unpack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf) = 0;
    virtual int getExtent() = 0;
    virtual int getSize() = 0;
    virtual void print(std::string indent) = 0;

    void (*packer)(void*, int, void*);
    void (*unpacker)(void*, int, void*);

    llvm::Function* FPack;
    llvm::Function* FUnpack;
};

/* Class for primitive types, such as MPI_INT, MPI_BYTE, etc */
class PrimitiveDatatype : public Datatype {
    public:
    enum PrimitiveType { BYTE, CHAR, DOUBLE, INT };   

    private:
    PrimitiveDatatype::PrimitiveType Type;
    int Extent;
    int Size;

    public:
    PrimitiveDatatype(PrimitiveType type);
    void Codegen_Pack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void Codegen_Unpack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    int getExtent();
    int getSize();
    void print(std::string indent);

};

/* Class for contiguous types */
class ContiguousDatatype : public Datatype {

    Datatype* Basetype;
    int Count;
    void Codegen(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf, int elemstride_in, int elemstride_out, bool pack);

    public:
    ContiguousDatatype(Datatype* type, int count) : Datatype(), Basetype(type), Count(count) {}
    void Codegen_Pack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void Codegen_Unpack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    int getExtent();
    int getSize();
    void print(std::string indent);

};

/* Class for vector types */
class VectorDatatype : public Datatype {

    Datatype* Basetype;
    int Count;
    int Blocklen;
    int Stride;

    public:
    VectorDatatype(Datatype* type, int count, int blocklen, int stride) : Datatype(), Basetype(type), Count(count), Blocklen(blocklen), Stride(stride) {}
    void Codegen_Pack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void Codegen_Unpack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    int getExtent();
    int getSize();
    void print(std::string indent);

};

/* Class for hvector types */
class HVectorDatatype : public Datatype {

    Datatype* Basetype;
    int Count;
    int Blocklen;
    int Stride;

    public:
    HVectorDatatype(Datatype* type, int count, int blocklen, int stride) :  Datatype(), Basetype(type), Count(count), Blocklen(blocklen), Stride(stride) {}
    void Codegen_Pack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void Codegen_Unpack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    int getExtent();
    int getSize();
    void print(std::string indent);

};

/* Class for indexed block types */
class IndexedBlockDatatype : public Datatype {

    int Count;
    Datatype* Basetype;
    int Blocklen;
    std::vector<int> Displ;
    void Codegen(llvm::Value *contig_buf, llvm::Value *hindexed_buf, llvm::Value* incount, bool gather);

    public:
    IndexedBlockDatatype(int count, int blocklen, int* displ, Datatype* basetype);
    void Codegen_Pack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void Codegen_Unpack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    int getExtent();
    int getSize();
    void print(std::string indent);

};

/* Class for hindexed types */
class HIndexedDatatype : public Datatype {

    int Count;
    Datatype* Basetype;
    std::vector<int> Blocklen;
    std::vector<int> Displ;
    void Codegen(llvm::Value *contig_buf, llvm::Value *hindexed_buf, llvm::Value* incount, bool gather);

    public:
    HIndexedDatatype(int count, int* blocklen, long* displ, Datatype* basetype);
    void Codegen_Pack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void Codegen_Unpack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    int getExtent();
    int getSize();
    void print(std::string indent);

};

/* Class for struct types */
class StructDatatype : public Datatype {

    int Count;
    std::vector<Datatype*> Types;
    std::vector<int> Blocklen;
    std::vector<int> Displ;

    public:
    StructDatatype(int count, int* blocklen, long*  displ, Datatype** types);
    void Codegen_Pack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void Codegen_Unpack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void Codegen(llvm::Value *compactbuf, llvm::Value *scatteredbuf, llvm::Value* incount, bool pack);
    int getExtent();
    int getSize();
    void print(std::string indent);

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

