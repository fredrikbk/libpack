#include <map>
#include <queue>
#include <cstdio>

#include "ddt_jit.hpp"

extern "C" {
#include "pack.h"
}

int LPK_Init() {
    farc::DDT_Init();
    return 0;
}

int LPK_Finalize() {
    farc::DDT_Finalize();
    return 0;
}

int LPK_Primitive(LPK_Primitivetype primitivetype, LPK_Datatype *newtype_p) {

    farc::Datatype* ddt = new farc::PrimitiveDatatype(farc::PrimitiveDatatype::INT);

    switch (primitivetype) {
        case LPK_CHAR:
            ddt = new farc::PrimitiveDatatype(farc::PrimitiveDatatype::CHAR);
            break;
        case LPK_BYTE:
            ddt = new farc::PrimitiveDatatype(farc::PrimitiveDatatype::BYTE);
            break;
        case LPK_DOUBLE:
            ddt = new farc::PrimitiveDatatype(farc::PrimitiveDatatype::DOUBLE);
            break;
        case LPK_FLOAT:
            ddt = new farc::PrimitiveDatatype(farc::PrimitiveDatatype::FLOAT);
            break;
        case LPK_INT:
            ddt = new farc::PrimitiveDatatype(farc::PrimitiveDatatype::INT);
            break;
        default:
            fprintf(stderr, "Error in %s line %i: primitive datatype not implemented\n", __FILE__, __LINE__);
    }

    *newtype_p = (void*) ddt;

    return 0;

}

int LPK_Contiguous(int count, LPK_Datatype oldtype, LPK_Datatype *newtype_p) {

    farc::Datatype* ntype = new farc::ContiguousDatatype(count, reinterpret_cast<farc::Datatype*>(oldtype));
    *newtype_p = reinterpret_cast<LPK_Datatype>(ntype);

    return 0;

}

int LPK_Vector(int count, int blocklen, int stride, LPK_Datatype oldtype, LPK_Datatype *newtype_p) {

    farc::Datatype* ntype = new farc::VectorDatatype(count, blocklen, stride,
                                 reinterpret_cast<farc::Datatype*>(oldtype)); 
    *newtype_p = reinterpret_cast<LPK_Datatype>(ntype);

    return 0;    

}

int LPK_Hvector(int count, int blocklen, LPK_Aint stride, LPK_Datatype oldtype, LPK_Datatype *newtype_p) {

    farc::Datatype* ntype = new farc::HVectorDatatype(count, blocklen, stride,
                                 reinterpret_cast<farc::Datatype*>(oldtype)); 
    *newtype_p = reinterpret_cast<LPK_Datatype>(ntype);

    return 0;
}

int LPK_Indexed_block(int count, int blocklen, int displacements[], LPK_Datatype oldtype, LPK_Datatype *newtype_p) {

    farc::Datatype* ntype = new farc::IndexedBlockDatatype(count, blocklen, displacements,
                                 reinterpret_cast<farc::Datatype*>(oldtype)); 
    *newtype_p = reinterpret_cast<LPK_Datatype>(ntype);

    return 0;

}

int LPK_Hindexed(int count, int blocklens[], LPK_Aint displacements[], LPK_Datatype oldtype, LPK_Datatype *newtype_p) {

    farc::Datatype* ntype = new farc::HIndexedDatatype(count, blocklens, displacements,
                                 reinterpret_cast<farc::Datatype*>(oldtype)); 
    *newtype_p = reinterpret_cast<LPK_Datatype>(ntype);

    return 0;

}

int LPK_Struct(int count, int blocklens[], LPK_Aint displacements[], LPK_Datatype oldtypes[], LPK_Datatype *newtype_p) {

    farc::Datatype* ntype = new farc::StructDatatype(count, blocklens, displacements,
                                 reinterpret_cast<farc::Datatype**>(oldtypes)); 
    *newtype_p = reinterpret_cast<LPK_Datatype>(ntype);

    return 0;

}

int LPK_Free(LPK_Datatype *ddt) {

    farc::DDT_Free(reinterpret_cast<farc::Datatype*>(ddt));

    return 0;

}

int LPK_Compile(LPK_Datatype *ddt) {

    farc::DDT_Commit(reinterpret_cast<farc::Datatype*>(ddt)); 

    return 0;

}

int LPK_Pack(void* inbuf, int incount, LPK_Datatype intype, void* outbuf) {

    farc::DDT_Unpack(inbuf, outbuf, reinterpret_cast<farc::Datatype*>(intype), incount);

    return 0;

}

int LPK_Unpack(void* inbuf, void* outbuf, int outcount, LPK_Datatype outtype) {

    farc::DDT_Unpack(inbuf, outbuf, reinterpret_cast<farc::Datatype*>(outtype), outcount);

    return 0;

}

int LPK_Get_extent(LPK_Datatype datatype, LPK_Aint *lb, LPK_Aint *extent) {

    *extent = reinterpret_cast<farc::Datatype*>(datatype)->getExtent();
    *lb = reinterpret_cast<farc::Datatype*>(datatype)->getLowerBound();

    return 0; 

}

int LPK_Get_size(LPK_Datatype datatype, int *size) {

    *size = reinterpret_cast<farc::Datatype*>(datatype)->getSize();

    return 0;
}



