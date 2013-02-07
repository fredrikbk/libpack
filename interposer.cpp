/**
The NAS DDT patches call the following functions:

MPI_Type_commit
MPI_Type_create_hindexed
MPI_Type_create_hvector
MPI_Type_create_struct
MPI_Type_f2c
MPI_Type_free
MPI_Type_vector

Further we need:
MPI_Irecv
MPI_Isend
MPI_Recv
MPI_Send
*/

#include <mpi.h>
#include "ddt_jit.hpp"

extern "C" {
int MPI_Type_commit(MPI_Datatype *datatype);
int MPI_Type_vector(int count, int blocklength, int stride, MPI_Datatype oldtype, MPI_Datatype *newtype);
int MPI_Type_hvector(int count, int blocklength, MPI_Aint stride, MPI_Datatype oldtype, MPI_Datatype *newtype);
int MPI_Type_hindexed(int count, int *array_of_blocklengths, MPI_Aint *array_of_displacements, MPI_Datatype oldtype, MPI_Datatype *newtype);
int MPI_Type_create_struct(int count, int array_of_blocklengths[], MPI_Aint array_of_displacements[], MPI_Datatype array_of_types[], MPI_Datatype *newtype);
int MPI_Type_free(MPI_Datatype *datatype);
//MPI_Datatype MPI_Type_f2c(MPI_Fint datatype);
int MPI_Isend(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request);
int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request);
int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status);
int MPI_Send(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
int MPI_Init(int *argc, char ***argv);
}

// we rely on MPI_Datatype beeing typedef'd to int so that we can use them as index into a vector
std::vector<FARC_Datatype> g_my_types;

int MPI_Init(int *argc, char ***argv) {

    /* we handle MPI_DOUBLE as a special case here, it's the only
       primitive datatype used in the NAS benchmarks, so MPI_DOUBLE
       is index 0 in our translation table - maybe we really should
       have the primitive ddts built in
     */

    FARC_PrimitiveDatatype* ddt = new  FARC_PrimitiveDatatype(MPI_DOUBLE);
    g_my_types.push_back(ddt);

    return PMPI_Init(argc, argv);

}

int MPI_Type_vector(int count, int blocklength, int stride, MPI_Datatype oldtype, MPI_Datatype *newtype) {

    FARC_Datatype* oldtype_farc;

    if (oldtype == MPI_DOUBLE) oldtype_farc = g_my_types[0];
    else oldtype_farc = g_my_types[oldtype];
    FARC_PrimitiveDatatype* ddt = new FARC_VectorDatatype(oldtype_farc, count, blocklength, stride);
    g_my_types.push_back(ddt);
    int idx = g_my_types.size() - 1;
    *newtype = idx;

    return MPI_SUCCESS;
}


int MPI_Type_commit(MPI_Datatype *datatype) {

    FARC_DDT_Commit(g_my_types[*((int*)(datatype))]);

    return MPI_SUCCESS;

}

int MPI_Send(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm) {

    if (datatype != MPI_DOUBLE) {
        int outsize = datatype->getSize() * count;
        char* outbuf = malloc(outsize);

        int (*FP)(void*, int, void*) = (int (*)(void*, int, void*))(intptr_t) farc_ddt->packer;
        FP(buf, count, outbuf);
    
        PMPI_Send(outbuf, outsize, MPI_BYTE, dest, tag, comm);

        free(outbuf);

        return MPI_SUCCESS;
    }
    else {
        PMPI_Send(outbuf, outsize, datatype, dest, tag, comm);
    }

}


