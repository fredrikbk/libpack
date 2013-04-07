#ifndef DDT_JIT_H
#define DDT_JIT_H

#include <cstdlib>
#include <vector>
#include <string>

#define LAZY           0
#define LLVM_OPTIMIZE  0 
#define VECTOR_UNROLL  0

#define TIME 0
#define DDT_OUTPUT     0 

// LLVM_OUTPUT should be picked up from the environment by the build system
#ifndef LLVM_OUTPUT
#define LLVM_OUTPUT 0
#endif

// If LLVM_OUTPUT is set then we always want lazy
#if LLVM_OUTPUT
#undef LAZY
#define LAZY 1
#endif
 

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
    virtual Datatype* Clone() = 0;
    virtual void Codegen_Pack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf) = 0;
    virtual void Codegen_Unpack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf) = 0;
    virtual int getExtent() = 0;
    virtual int getSize() = 0;
	virtual std::string toString() = 0;
    virtual void print();

    void (*packer)(void*, int, void*);
    void (*unpacker)(void*, int, void*);

    llvm::Function* FPack;
    llvm::Function* FUnpack;
};

/* Class for primitive types, such as MPI_INT, MPI_BYTE, etc */
class PrimitiveDatatype : public Datatype {
    public:
    enum PrimitiveType { BYTE, CHAR, DOUBLE, FLOAT, INT };   

    private:
    PrimitiveDatatype::PrimitiveType Type;
    int Extent;
    int Size;

    public:
    PrimitiveDatatype(PrimitiveType type);
    ~PrimitiveDatatype(void) {};
    void Codegen_Pack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void Codegen_Unpack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    int getExtent();
    int getSize();
    std::string toString();
    PrimitiveDatatype* Clone();
};

/* Class for contiguous types */
class ContiguousDatatype : public Datatype {

    Datatype* Basetype;
    int Count;
    void Codegen(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf, int elemstride_in, int elemstride_out, bool pack);

    public:
    ContiguousDatatype(Datatype* type, int count);
    ~ContiguousDatatype(void);
    void Codegen_Pack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void Codegen_Unpack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    int getExtent();
    int getSize();
    std::string toString();
    ContiguousDatatype* Clone();

};

/* Class for vector types */
class VectorDatatype : public Datatype {

    Datatype* Basetype;
    int Count;
    int Blocklen;
    int Stride;

    public:
    VectorDatatype(Datatype* type, int count, int blocklen, int stride); 
    ~VectorDatatype(void); 
    void Codegen_Pack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void Codegen_Unpack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    int getExtent();
    int getSize();
    std::string toString();
    VectorDatatype* Clone();

};

/* Class for hvector types */
class HVectorDatatype : public Datatype {

    Datatype* Basetype;
    int Count;
    int Blocklen;
    int Stride;

    public:
    HVectorDatatype(Datatype* type, int count, int blocklen, int stride);
    ~HVectorDatatype(void);
    void Codegen_Pack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void Codegen_Unpack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    int getExtent();
    int getSize();
    std::string toString();
    HVectorDatatype* Clone();

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
    ~IndexedBlockDatatype(void);
    void Codegen_Pack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void Codegen_Unpack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    int getExtent();
    int getSize();
    std::string toString();
    IndexedBlockDatatype* Clone();

};

/* Class for hindexed types */
class HIndexedDatatype : public Datatype {

    int Count;
    Datatype* Basetype;
    std::vector<int> Blocklen;
    std::vector<long> Displ;
    void Codegen(llvm::Value *contig_buf, llvm::Value *hindexed_buf, llvm::Value* incount, bool gather);

    public:
    HIndexedDatatype(int count, int* blocklen, long* displ, Datatype* basetype);
    ~HIndexedDatatype(void);
    void Codegen_Pack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void Codegen_Unpack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    int getExtent();
    int getSize();
    std::string toString();
    HIndexedDatatype* Clone();

};

/* Class for struct types */
class StructDatatype : public Datatype {

    int Count;
    std::vector<Datatype*> Types;
    std::vector<int> Blocklen;
    std::vector<long> Displ;

    public:
    StructDatatype(int count, int* blocklen, long*  displ, Datatype** types);
    ~StructDatatype(void);
    void Codegen_Pack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void Codegen_Unpack(llvm::Value* inbuf, llvm::Value* incount, llvm::Value* outbuf);
    void Codegen(llvm::Value *compactbuf, llvm::Value *scatteredbuf, llvm::Value* incount, bool pack);
    int getExtent();
    int getSize();
    std::string toString();
    StructDatatype* Clone();

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

