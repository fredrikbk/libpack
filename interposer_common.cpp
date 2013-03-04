#include "interposer_common.h"
#include <stdio.h>
#include "ddt_jit.hpp"


struct recvop {
    void* buf;
    void* shaddow;
    int count;
    MPI_Datatype datatype;
};

static std::map<MPI_Datatype, FARC_Datatype*> g_my_types;
static std::map<MPI_Request, char*> g_my_buffers;
static std::map<MPI_Request, struct recvop> g_my_recv_reqs;

void interposer_init() {
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
}

void interposer_hvector(int count, int blocklength, MPI_Aint stride, MPI_Datatype oldtype, MPI_Datatype *newtype) {
    FARC_Datatype* oldtype_farc;
    oldtype_farc = g_my_types[oldtype];
    FARC_Datatype* ddt = new FARC_HVectorDatatype(oldtype_farc, count, blocklength, stride);
    g_my_types[*newtype] = ddt;
}

void interposer_vector(int count, int blocklength, int stride, MPI_Datatype oldtype, MPI_Datatype *newtype) {
    FARC_Datatype* oldtype_farc;
    oldtype_farc = g_my_types[oldtype];
    FARC_Datatype* ddt = new FARC_VectorDatatype(oldtype_farc, count, blocklength, stride);
    g_my_types[*newtype] = ddt;
}

void interposer_create_struct(int count, int array_of_blocklengths[], MPI_Aint array_of_displacements[], MPI_Datatype array_of_types[], MPI_Datatype *newtype) {
    FARC_Datatype** farc_oldtypes = (FARC_Datatype**) malloc(count * sizeof(FARC_Datatype*));
    for (int i=0; i<count; i++) {
        farc_oldtypes[i] = g_my_types[array_of_types[i]];
    }
    FARC_Datatype* ddt = new FARC_StructDatatype(count, array_of_blocklengths, array_of_displacements, farc_oldtypes);
    g_my_types[*newtype] = ddt;
}

void interposer_struct(int count, int *array_of_blocklengths, MPI_Aint *array_of_displacements, MPI_Datatype *array_of_types, MPI_Datatype *newtype) {
    FARC_Datatype** farc_oldtypes = (FARC_Datatype**) malloc(count * sizeof(FARC_Datatype*));
    for (int i=0; i<count; i++) {
        farc_oldtypes[i] = g_my_types[array_of_types[i]];
    }
    FARC_Datatype* ddt = new FARC_StructDatatype(count, array_of_blocklengths, array_of_displacements, farc_oldtypes);
    free(farc_oldtypes);
    g_my_types[*newtype] = ddt;
}

void interposer_create_hvector(int count, int blocklength, MPI_Aint stride, MPI_Datatype oldtype, MPI_Datatype *newtype) {
    FARC_Datatype* oldtype_farc;
    oldtype_farc = g_my_types[oldtype];
    FARC_Datatype* ddt = new FARC_HVectorDatatype(oldtype_farc, count, blocklength, stride);
    g_my_types[*newtype] = ddt;
}

void interposer_create_hindexed(int count, int array_of_blocklengths[], MPI_Aint array_of_displacements[], MPI_Datatype oldtype, MPI_Datatype *newtype) {
    FARC_Datatype* oldtype_farc;
    oldtype_farc = g_my_types[oldtype];
    FARC_Datatype* ddt = new FARC_HIndexedDatatype(count, array_of_blocklengths, array_of_displacements, oldtype_farc);
    g_my_types[*newtype] = ddt;
}

void interposer_hindexed(int count, int array_of_blocklengths[], MPI_Aint array_of_displacements[], MPI_Datatype oldtype, MPI_Datatype *newtype) {
    FARC_Datatype* oldtype_farc;
    oldtype_farc = g_my_types[oldtype];
    FARC_Datatype* ddt = new FARC_HIndexedDatatype(count, array_of_blocklengths, array_of_displacements, oldtype_farc);
    g_my_types[*newtype] = ddt;
}

void interposer_contiguous(int count, MPI_Datatype oldtype, MPI_Datatype *newtype) {
    FARC_Datatype* oldtype_farc;
    oldtype_farc = g_my_types[oldtype];
    FARC_Datatype* ddt = new FARC_ContiguousDatatype(oldtype_farc, count);
    g_my_types[*newtype] = ddt;
}

void interposer_commit(MPI_Datatype *datatype) {
    FARC_DDT_Commit(g_my_types[*datatype]);
}

void interposer_free(MPI_Datatype *datatype) {
    FARC_DDT_Free(g_my_types[*datatype]);
}

// OPT: Allocate buffer on stack or cache buffers in a map keyed on sizes
char* interposer_buffer_alloc(int count, MPI_Datatype datatype, int* buf_size) {
    *buf_size = g_my_types[datatype]->getSize() * count;
    return (char*) malloc(*buf_size);
}

void interposer_buffer_free(char* buf) {
    free(buf);
}

void interposer_buffer_register(MPI_Request* request, char* buf) {
    g_my_buffers[*request] = buf;
}

void interposer_recvop_register(void *data, int count, MPI_Datatype datatype, char* buf, MPI_Request *request) {
    struct recvop rop;
    rop.buf = data;
    rop.count = count;
    rop.datatype = datatype;
    rop.shaddow = buf;
    g_my_recv_reqs[*request] = rop;
}

char* interposer_pack(void *data, int count, MPI_Datatype datatype, int *buf_size) {
    char* buf = interposer_buffer_alloc(count, datatype, buf_size);

    int (*FP)(void*, int, void*) = (int (*)(void*, int, void*))(intptr_t) g_my_types[datatype]->packer;
    FP(data, count, buf);

    return buf;
}

void interposer_unpack(void *data, int count, MPI_Datatype datatype, char* buf) {
    int (*FP)(void*, int, void*) = (int (*)(void*, int, void*))(intptr_t) g_my_types[datatype]->unpacker;
    FP(buf, count, data);
}

//**********************************************************

int MPI_Waitall(int count, MPI_Request *array_of_requests, MPI_Status *array_of_statuses) {

    // we need to copy the old requests here :-( (they will become MPI_REQUEST_NULL)
    // OPT: Allocate these on the stack
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





// Fortran Bindings
extern "C" {
    void interposer_init_() {
	interposer_init();
    }
    
    void interposer_hvector_(int count, int blocklength, MPI_Aint stride, MPI_Datatype oldtype, MPI_Datatype *newtype) {
	interposer_hvector(count, blocklength, stride, oldtype, newtype);
    }

    void interposer_vector_(int count, int blocklength, int stride, MPI_Datatype oldtype, MPI_Datatype *newtype) {
	interposer_vector(count, blocklength, stride, oldtype, newtype);
    }

    void interposer_create_struct_(int count, int array_of_blocklengths[], MPI_Aint array_of_displacements[], MPI_Datatype array_of_types[], MPI_Datatype *newtype) {
	interposer_create_struct(count, array_of_blocklengths, array_of_displacements, array_of_types, newtype);
    }
	
    void interposer_struct_(int count, int *array_of_blocklengths, MPI_Aint *array_of_displacements, MPI_Datatype *array_of_types, MPI_Datatype *newtype) {
	interposer_struct(count, array_of_blocklengths, array_of_displacements, array_of_types, newtype);
    }

    void interposer_create_hvector_(int count, int blocklength, MPI_Aint stride, MPI_Datatype oldtype, MPI_Datatype *newtype) {
	interposer_create_hvector(count, blocklength, stride, oldtype, newtype);
    }

    void interposer_create_hindexed_(int count, int array_of_blocklengths[], MPI_Aint array_of_displacements[], MPI_Datatype oldtype, MPI_Datatype *newtype) {
	interposer_create_hindexed(count, array_of_blocklengths, array_of_displacements, oldtype, newtype);
    }

    void interposer_hindexed_(int count, int array_of_blocklengths[], MPI_Aint array_of_displacements[], MPI_Datatype oldtype, MPI_Datatype *newtype) {
	interposer_hindexed(count, array_of_blocklengths, array_of_displacements, oldtype, newtype);
    }

    void interposer_contiguous_(int count, MPI_Datatype oldtype, MPI_Datatype *newtype) {
	interposer_contiguous(count, oldtype, newtype);
    }

    void interposer_commit_(MPI_Datatype *datatype) {
	interposer_commit(datatype);
    }

    void interposer_free_(MPI_Datatype *datatype) {
	interposer_free(datatype);
    }

    char* interposer_buffer_alloc_(int count, MPI_Datatype datatype, int* buf_size) {
	return interposer_buffer_alloc(count, datatype, buf_size);
    }

    void interposer_buffer_free_(char* buf) {
	interposer_buffer_free(buf);
    }

    void interposer_buffer_register_(MPI_Request* request, char* buf) {
	interposer_buffer_register(request, buf);
    }

    void interposer_recvop_register_(void *data, int count, MPI_Datatype datatype, char* buf, MPI_Request *request) {
	interposer_recvop_register(data, count, datatype, buf, request);
    }

    char* interposer_pack_(void *data, int count, MPI_Datatype datatype, int *buf_size) {
	return interposer_pack(data, count, datatype, buf_size);
    }

    void interposer_unpack_(void *data, int count, MPI_Datatype datatype, char* buf) {
	interposer_unpack(data, count, datatype, buf);
    }
}
