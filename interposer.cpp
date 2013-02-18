/**
Link this into your application to use FARC DDTs (see the pmpi-tests Makefile for how to do that)
*/

#include <mpi.h>
#include "ddt_jit.hpp"

struct recvop {
    void* buf;
    void* shaddow;
    int count;
    MPI_Datatype datatype;
};

std::map<MPI_Datatype, FARC_Datatype*> g_my_types;
std::map<MPI_Request, char*> g_my_buffers;
std::map<MPI_Request, struct recvop> g_my_recv_reqs;

int MPI_Init(int *argc, char ***argv) {

    FARC_PrimitiveDatatype* ddt;

    ddt = new FARC_PrimitiveDatatype(MPI_DOUBLE);
    g_my_types[MPI_DOUBLE] = ddt;

    ddt = new FARC_PrimitiveDatatype(MPI_INT);
    g_my_types[MPI_INT] = ddt;

    ddt = new FARC_PrimitiveDatatype(MPI_BYTE);
    g_my_types[MPI_BYTE] = ddt;

    ddt = new FARC_PrimitiveDatatype(MPI_CHAR);
    g_my_types[MPI_CHAR] = ddt;

    //TODO add other primitive types here

    FARC_DDT_Init();

    return PMPI_Init(argc, argv);

}

int MPI_Type_hvector(int count, int blocklength, MPI_Aint stride, MPI_Datatype oldtype, MPI_Datatype *newtype) {

    PMPI_Type_hvector(count, blocklength, stride, oldtype, newtype);

    FARC_Datatype* oldtype_farc;
    oldtype_farc = g_my_types[oldtype];
    FARC_Datatype* ddt = new FARC_HVectorDatatype(oldtype_farc, count, blocklength, stride);

    g_my_types[*newtype] = ddt;

    return MPI_SUCCESS;

}

int MPI_Type_vector(int count, int blocklength, int stride, MPI_Datatype oldtype, MPI_Datatype *newtype) {

    PMPI_Type_vector(count, blocklength, stride, oldtype, newtype);

    FARC_Datatype* oldtype_farc;
    oldtype_farc = g_my_types[oldtype];
    FARC_Datatype* ddt = new FARC_VectorDatatype(oldtype_farc, count, blocklength, stride);

    g_my_types[*newtype] = ddt;

    return MPI_SUCCESS;

}

int MPI_Type_create_struct(int count, int array_of_blocklengths[], MPI_Aint array_of_displacements[], MPI_Datatype array_of_types[], MPI_Datatype *newtype) {

    PMPI_Type_create_struct(count, array_of_blocklengths, array_of_displacements, array_of_types, newtype);

    FARC_Datatype** farc_oldtypes = (FARC_Datatype**) malloc(count * sizeof(FARC_Datatype*));
    for (int i=0; i<count; i++) {
        farc_oldtypes[i] = g_my_types[array_of_types[i]];
    }
    FARC_Datatype* ddt = new FARC_StructDatatype(count, array_of_blocklengths, array_of_displacements, farc_oldtypes);
    g_my_types[*newtype] = ddt;

    return MPI_SUCCESS;

}

int MPI_Type_struct(int count, int *array_of_blocklengths, MPI_Aint *array_of_displacements, MPI_Datatype *array_of_types, MPI_Datatype *newtype) {

    PMPI_Type_struct(count, array_of_blocklengths, array_of_displacements, array_of_types, newtype);

    FARC_Datatype** farc_oldtypes = (FARC_Datatype**) malloc(count * sizeof(FARC_Datatype*));
    for (int i=0; i<count; i++) {
        farc_oldtypes[i] = g_my_types[array_of_types[i]];
    }
    FARC_Datatype* ddt = new FARC_StructDatatype(count, array_of_blocklengths, array_of_displacements, farc_oldtypes);
    free(farc_oldtypes);
    g_my_types[*newtype] = ddt;

    return MPI_SUCCESS;

}

int MPI_Type_create_hvector(int count, int blocklength, MPI_Aint stride, MPI_Datatype oldtype, MPI_Datatype *newtype) {

    PMPI_Type_create_hvector(count, blocklength, stride, oldtype, newtype);

    FARC_Datatype* oldtype_farc;
    oldtype_farc = g_my_types[oldtype];
    FARC_Datatype* ddt = new FARC_HVectorDatatype(oldtype_farc, count, blocklength, stride);

    g_my_types[*newtype] = ddt;

    return MPI_SUCCESS;

}

int MPI_Type_create_hindexed(int count, int array_of_blocklengths[], MPI_Aint array_of_displacements[], MPI_Datatype oldtype, MPI_Datatype *newtype) {

    PMPI_Type_create_hindexed(count, array_of_blocklengths, array_of_displacements, oldtype, newtype);

    FARC_Datatype* oldtype_farc;
    oldtype_farc = g_my_types[oldtype];
    FARC_Datatype* ddt = new FARC_HIndexedDatatype(count, array_of_blocklengths, array_of_displacements, oldtype_farc);

    g_my_types[*newtype] = ddt;

    return MPI_SUCCESS;

}

int MPI_Type_hindexed(int count, int array_of_blocklengths[], MPI_Aint array_of_displacements[], MPI_Datatype oldtype, MPI_Datatype *newtype) {

    PMPI_Type_hindexed(count, array_of_blocklengths, array_of_displacements, oldtype, newtype);

    FARC_Datatype* oldtype_farc;
    oldtype_farc = g_my_types[oldtype];
    FARC_Datatype* ddt = new FARC_HIndexedDatatype(count, array_of_blocklengths, array_of_displacements, oldtype_farc);

    g_my_types[*newtype] = ddt;

    return MPI_SUCCESS;

}

int MPI_Type_contiguous(int count, MPI_Datatype oldtype, MPI_Datatype *newtype) {

    PMPI_Type_contiguous(count, oldtype, newtype);

    FARC_Datatype* oldtype_farc;
    oldtype_farc = g_my_types[oldtype];
    FARC_Datatype* ddt = new FARC_ContiguousDatatype(oldtype_farc, count);

    g_my_types[*newtype] = ddt;

    return MPI_SUCCESS;

}

int MPI_Type_commit(MPI_Datatype *datatype) {

    PMPI_Type_commit(datatype);
    FARC_DDT_Commit(g_my_types[*datatype]);

    return MPI_SUCCESS;

}

int MPI_Type_free(MPI_Datatype *datatype) {

    FARC_DDT_Free(g_my_types[*datatype]);
    int ret = PMPI_Type_free(datatype);

    return ret;
}

int MPI_Send(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm) {

    // TODO check other primitive types!
    if (datatype != MPI_DOUBLE) {
        int outsize = g_my_types[datatype]->getSize() * count;
        char* outbuf = (char*) malloc(outsize);

        int (*FP)(void*, int, void*) = (int (*)(void*, int, void*))(intptr_t) g_my_types[datatype]->packer;
        FP(buf, count, outbuf);
    
        PMPI_Send(outbuf, outsize, MPI_BYTE, dest, tag, comm);

        free(outbuf);

        return MPI_SUCCESS;
    }
    else {
        PMPI_Send(buf, count, MPI_BYTE, dest, tag, comm);
    }

    return MPI_SUCCESS;
}

int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status) {

    // TODO check other primitive types!


    if (datatype != MPI_DOUBLE) {
        int insize = g_my_types[datatype]->getSize() * count;
        char* inbuf = (char*) malloc(insize);

        PMPI_Recv(inbuf, insize, MPI_BYTE, source, tag, comm, status);
        int (*FP)(void*, int, void*) = (int (*)(void*, int, void*))(intptr_t) g_my_types[datatype]->unpacker;
        FP(inbuf, count, buf);

        free(inbuf);

        return MPI_SUCCESS;
    }
    else {
        PMPI_Recv(buf, count, MPI_BYTE, source, tag, comm, status);
    }

    return MPI_SUCCESS;
   
}

int MPI_Isend(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request) {

    // TODO check other primitive types!
    if (datatype != MPI_DOUBLE) {
        int outsize = g_my_types[datatype]->getSize() * count;
        char* outbuf = (char*) malloc(outsize);

        int (*FP)(void*, int, void*) = (int (*)(void*, int, void*))(intptr_t) g_my_types[datatype]->packer;
        FP(buf, count, outbuf);
    
        PMPI_Isend(outbuf, outsize, MPI_BYTE, dest, tag, comm, request);

        g_my_buffers[*request] = outbuf;

        return MPI_SUCCESS;
    }
    else {
        PMPI_Isend(buf, count, MPI_BYTE, dest, tag, comm, request);
        g_my_buffers[*request] = NULL;
    }

    return MPI_SUCCESS;
}

int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request) {

    struct recvop rop;
    rop.buf = buf;
    rop.count = count;
    rop.datatype = datatype;

    // TODO check other primitive types!

    if (datatype != MPI_DOUBLE) {
        int insize = g_my_types[datatype]->getSize() * count;
        char* inbuf = (char*) malloc(insize);
        rop.shaddow = inbuf;

        PMPI_Irecv(inbuf, insize, MPI_BYTE, source, tag, comm, request);
        g_my_buffers[*request] = inbuf;
        g_my_recv_reqs[*request] = rop;

        return MPI_SUCCESS;
    }
    else {
        PMPI_Irecv(buf, count, MPI_BYTE, source, tag, comm, request);
        g_my_buffers[*request] = NULL;
    }

    return MPI_SUCCESS;
   
}

int MPI_Waitall(int count, MPI_Request *array_of_requests, MPI_Status *array_of_statuses) {

    // we need to copy the old requests here :-( (they will become MPI_REQUEST_NULL)
    MPI_Request* oldrequests = (MPI_Request*) malloc(count * sizeof(MPI_Request));
    for (int i=0; i<count; i++) {
        oldrequests[i] = array_of_requests[i];
    }

    int ret = PMPI_Waitall(count, array_of_requests, array_of_statuses);

    // now go over them and unpack recv requests and free buffers
    for (int i=0; i<count; i++) {
        if (oldrequests[i] != MPI_REQUEST_NULL) {
             //check if it was a recvop, if yes, unpack TODO this is expensive :-(
            std::map<MPI_Request, struct recvop>::iterator it;
            it = g_my_recv_reqs.find(oldrequests[i]);
            if (it != g_my_recv_reqs.end()) {
                int (*FP)(void*, int, void*) = (int (*)(void*, int, void*))(intptr_t) g_my_types[it->second.datatype]->unpacker;
                FP(it->second.shaddow, it->second.count, it->second.buf);
                g_my_recv_reqs.erase(it);
            }
            free(g_my_buffers[oldrequests[i]]);
            g_my_buffers.erase(oldrequests[i]);
        }
    }

    free(oldrequests);

    return ret;
}

int MPI_Testall(int count, MPI_Request *array_of_requests, int *flag, MPI_Status *array_of_statuses) {

    // we need to copy the old requests here :-( (they will become MPI_REQUEST_NULL)
    MPI_Request* oldrequests = (MPI_Request*) malloc(count * sizeof(MPI_Request));
    for (int i=0; i<count; i++) {
        oldrequests[i] = array_of_requests[i];
    }

    int ret = PMPI_Testall(count, array_of_requests, flag, array_of_statuses);

    // now go over them and unpack recv requests and free buffers
    for (int i=0; i<count; i++) {
        if ((oldrequests[i] != MPI_REQUEST_NULL) && (array_of_requests[i] == MPI_REQUEST_NULL)) {
             //check if it was a recvop, if yes, unpack TODO this is expensive :-(
            std::map<MPI_Request, struct recvop>::iterator it;
            it = g_my_recv_reqs.find(oldrequests[i]);
            if (it != g_my_recv_reqs.end()) {
                int (*FP)(void*, int, void*) = (int (*)(void*, int, void*))(intptr_t) g_my_types[it->second.datatype]->unpacker;
                FP(it->second.shaddow, it->second.count, it->second.buf);
                g_my_recv_reqs.erase(it);
            }
            free(g_my_buffers[oldrequests[i]]);
            g_my_buffers.erase(oldrequests[i]);
        }
    }

    free(oldrequests);

    return ret;

}

//TODO we should also implement testany, testsome, waitany, waitsome to be complete, but who uses those :D

int MPI_Wait(MPI_Request *request, MPI_Status *status) {

    int ret;

    if (*request != MPI_REQUEST_NULL) {
        MPI_Request oldrequest = *request;
        ret = PMPI_Wait(request, status);

        //check if it was a recvop, if yes, unpack TODO this is expensive :-(
        std::map<MPI_Request, struct recvop>::iterator it;
        it = g_my_recv_reqs.find(oldrequest);
        if (it != g_my_recv_reqs.end()) {
            int (*FP)(void*, int, void*) = (int (*)(void*, int, void*))(intptr_t) g_my_types[it->second.datatype]->unpacker;
            FP(it->second.shaddow, it->second.count, it->second.buf);
            g_my_recv_reqs.erase(it);
        }

        free(g_my_buffers[oldrequest]);
        g_my_buffers.erase(oldrequest);
    }
    else {
        ret = PMPI_Wait(request, status);
    }

    return ret;

}

int MPI_Test(MPI_Request *request, int *flag, MPI_Status *status) {

    int ret;

    if (*request != MPI_REQUEST_NULL) {
        MPI_Request oldrequest = *request;
        ret = PMPI_Test(request, flag, status);
        if (*request == MPI_REQUEST_NULL) {
             //check if it was a recvop, if yes, unpack TODO this is expensive :-(
            std::map<MPI_Request, struct recvop>::iterator it;
            it = g_my_recv_reqs.find(oldrequest);
            if (it != g_my_recv_reqs.end()) {
                int (*FP)(void*, int, void*) = (int (*)(void*, int, void*))(intptr_t) g_my_types[it->second.datatype]->unpacker;
                FP(it->second.shaddow, it->second.count, it->second.buf);
                g_my_recv_reqs.erase(it);
            }
           
            free(g_my_buffers[oldrequest]);
            g_my_buffers.erase(oldrequest);
        }
    }
    else {
        ret = PMPI_Test(request, flag, status);
    }

    return ret;

}




