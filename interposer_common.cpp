#include "interposer_common.h"

#include <cstdlib>
#include <string>
#include <map>
#include <queue>
#include <list>
#include <mpi.h>

#include "ddt_jit.hpp"

//#include "copy_benchmark/hrtimer/hrtimer.h"
//static HRT_TIMESTAMP_T start, stop;
//static uint64_t tmp;

using namespace farc;

struct Request {
    MPI_Request *mpi_req;

    void* tmpbuf;

    // If it is a Irecv op
    void* usrbuf;
    int count;
    MPI_Datatype datatype;
};

static inline bool is_recv(Request req) {
    return req.usrbuf != NULL;
}

static std::list<struct Request> g_outstanding_requests; 

/* Datatype lookup data structures */
#define DDT_FAST_CACHE_SIZE 50
static Datatype* g_types[DDT_FAST_CACHE_SIZE];
static std::queue<int> g_types_freelist;
static std::map<MPI_Datatype, Datatype*> g_types_fallback;

static PrimitiveDatatype farc_double(PrimitiveDatatype::DOUBLE);
static PrimitiveDatatype farc_float(PrimitiveDatatype::FLOAT);
static PrimitiveDatatype farc_int(PrimitiveDatatype::INT);
static PrimitiveDatatype farc_byte(PrimitiveDatatype::BYTE);
static PrimitiveDatatype farc_char(PrimitiveDatatype::CHAR);
//TODO add other primitive types here

static inline MPI_Datatype datatype_handle_create() {
    MPI_Datatype ddt = (MPI_Datatype)g_types_freelist.front();
    g_types_freelist.pop();
    return ddt;
}

static inline void datatype_handle_free(MPI_Datatype* ddt_handle) {
    g_types_freelist.push(*ddt_handle);
}

static inline void datatype_store(MPI_Datatype dt_handle, Datatype *dt) {
    if ((int)dt_handle < 50) {
        g_types[(int)dt_handle] = dt;
    }
    else {
        g_types_fallback[dt_handle] = dt;
    }
}

static inline Datatype* datatype_retrieve(MPI_Datatype dt_handle) {

    switch (dt_handle) {
        case MPI_DOUBLE:
            return &farc_double;
            break;
        case MPI_INT:
            return &farc_int;
            break;
        case MPI_BYTE:
            return &farc_byte;
            break;
        case MPI_CHAR:
            return &farc_char;
            break;
        case MPI_FLOAT:
            return &farc_float;
            break;
    }

    return ((int)dt_handle < 50) ? g_types[(int)dt_handle]
                                  : g_types_fallback[dt_handle];
}

void interposer_init() {
    // Populate the free list (primitive datatype locations become dead slots)
    // Note: MPICH has no primitive datatypes in the range 0:49, but that might
    //       not hold for other MPI implementations or the future so a more 
    //       robust solution should remove all primitive datatypes from the 
    //       freelist if they are there
    for (int i=0; i < DDT_FAST_CACHE_SIZE; ++i) {
        g_types_freelist.push(i);
    } 

    DDT_Init();
}

void interposer_finalize() {
    DDT_Finalize();
}

void interposer_hvector(int count, int blocklength, MPI_Aint stride, MPI_Datatype oldtype, MPI_Datatype *newtype) {

    *newtype = datatype_handle_create();

    Datatype* oldtype_farc;
    oldtype_farc = datatype_retrieve(oldtype);
    Datatype* ddt = new HVectorDatatype(oldtype_farc, count, blocklength, stride);
    datatype_store(*newtype, ddt);
}

void interposer_vector(int count, int blocklength, int stride, MPI_Datatype oldtype, MPI_Datatype *newtype) {

    *newtype = datatype_handle_create();

    Datatype* oldtype_farc;
    oldtype_farc = datatype_retrieve(oldtype);
    Datatype* ddt = new VectorDatatype(oldtype_farc, count, blocklength, stride);
    datatype_store(*newtype, ddt);
}

void interposer_create_struct(int count, int array_of_blocklengths[], MPI_Aint array_of_displacements[], MPI_Datatype array_of_types[], MPI_Datatype *newtype) {

    *newtype = datatype_handle_create();

    Datatype** farc_oldtypes = (Datatype**) malloc(count * sizeof(Datatype*));
    for (int i=0; i<count; i++) {
        farc_oldtypes[i] = datatype_retrieve(array_of_types[i]);
    }
    Datatype* ddt = new StructDatatype(count, array_of_blocklengths, array_of_displacements, farc_oldtypes);
    datatype_store(*newtype, ddt);
}

void interposer_struct(int count, int *array_of_blocklengths, MPI_Aint *array_of_displacements, MPI_Datatype *array_of_types, MPI_Datatype *newtype) {

    *newtype = datatype_handle_create();

    Datatype** farc_oldtypes = (Datatype**) malloc(count * sizeof(Datatype*));
    for (int i=0; i<count; i++) {
        farc_oldtypes[i] = datatype_retrieve(array_of_types[i]);
    }
    Datatype* ddt = new StructDatatype(count, array_of_blocklengths, array_of_displacements, farc_oldtypes);
    free(farc_oldtypes);
    datatype_store(*newtype, ddt);
}

void interposer_create_hvector(int count, int blocklength, MPI_Aint stride, MPI_Datatype oldtype, MPI_Datatype *newtype) {

    *newtype = datatype_handle_create();

    Datatype* oldtype_farc;
    oldtype_farc = datatype_retrieve(oldtype);
    Datatype* ddt = new HVectorDatatype(oldtype_farc, count, blocklength, stride);
    datatype_store(*newtype, ddt);
}

void interposer_create_hindexed(int count, int array_of_blocklengths[], MPI_Aint array_of_displacements[], MPI_Datatype oldtype, MPI_Datatype *newtype) {

    *newtype = datatype_handle_create();

    Datatype* oldtype_farc;
    oldtype_farc = datatype_retrieve(oldtype);
    Datatype* ddt = new HIndexedDatatype(count, array_of_blocklengths, array_of_displacements, oldtype_farc);
    datatype_store(*newtype, ddt);
}

void interposer_hindexed(int count, int array_of_blocklengths[], MPI_Aint array_of_displacements[], MPI_Datatype oldtype, MPI_Datatype *newtype) {

    *newtype = datatype_handle_create();

    Datatype* oldtype_farc;
    oldtype_farc = datatype_retrieve(oldtype);
    Datatype* ddt = new HIndexedDatatype(count, array_of_blocklengths, array_of_displacements, oldtype_farc);
    datatype_store(*newtype, ddt);
}

void interposer_indexed_block(int count, int blocklength, int array_of_displacements[], MPI_Datatype oldtype, MPI_Datatype *newtype) {

    *newtype = datatype_handle_create();

    Datatype* oldtype_farc;
    oldtype_farc = datatype_retrieve(oldtype);
    Datatype* ddt = new IndexedBlockDatatype(count, blocklength, array_of_displacements, oldtype_farc);
    datatype_store(*newtype, ddt);
    
}

void interposer_contiguous(int count, MPI_Datatype oldtype, MPI_Datatype *newtype) {

    *newtype = datatype_handle_create();

    Datatype* oldtype_farc;
    oldtype_farc = datatype_retrieve(oldtype);
    Datatype* ddt = new ContiguousDatatype(oldtype_farc, count);
    datatype_store(*newtype, ddt);
}

void interposer_commit(MPI_Datatype *datatype) {
    DDT_Commit(datatype_retrieve(*datatype));
}

void interposer_free(MPI_Datatype *datatype) {
    DDT_Free(datatype_retrieve(*datatype));
    datatype_handle_free(datatype);
}

int interposer_type_size(MPI_Datatype datatype) {
    return datatype_retrieve(datatype)->getSize();
}

int interposer_type_extent(MPI_Datatype datatype) {
    return datatype_retrieve(datatype)->getExtent();
}


// TODO: Implement segmenting (and possibly reduce buffer size)
const int scratch_size = 2 * 1024 * 1024;
bool scratch_in_use = false;
static char scratch[scratch_size];

void* interposer_buffer_alloc(int count, MPI_Datatype datatype, int* buf_size) {
    *buf_size = datatype_retrieve(datatype)->getSize() * count;

    if (*buf_size <= scratch_size && !scratch_in_use) {
        scratch_in_use = true;
        return scratch;
    }
    else {
        return malloc(*buf_size);
    }
}

void interposer_buffer_free(void* tmpbuf) {
    if (tmpbuf != scratch) {
        free(tmpbuf);
    }
    else {
        scratch_in_use = true;
    }
}

void* interposer_pack(void *data, int count, MPI_Datatype datatype, int *buf_size) {
    void* buf = interposer_buffer_alloc(count, datatype, buf_size);
    DDT_Pack(data, buf, datatype_retrieve(datatype), count);
    return buf;
}

void interposer_pack_providedbuf(void* inbuf, int incount, MPI_Datatype datatype, void *outbuf) {
    DDT_Pack(inbuf, outbuf, datatype_retrieve(datatype), incount);
}

void interposer_unpack(void *data, int count, MPI_Datatype datatype, void* buf) {
    DDT_Unpack(buf, data, datatype_retrieve(datatype), count);
}

void interposer_request_register(void *tmpbuf, void *usrbuf, int count, MPI_Datatype datatype, MPI_Request *request) {
    struct Request req;
    req.mpi_req = request;

    req.tmpbuf = tmpbuf;

    req.usrbuf = usrbuf;
    req.count = count;
    req.datatype = datatype;
    g_outstanding_requests.push_back(req);

    if (tmpbuf != NULL && is_recv(req)) {
        DDT_Lazy_Unpack_Commit(datatype_retrieve(datatype));
    }
}

void interposer_request_free(MPI_Request *request) {
    for (std::list<Request>::iterator req = g_outstanding_requests.begin();
            req != g_outstanding_requests.end(); req++) {
        if (req->mpi_req == request) {
            if (req->tmpbuf != NULL) {
                // If it was a recv request then unpack it 
                if (is_recv(*req)) {
                    DDT_Unpack(req->tmpbuf, req->usrbuf, datatype_retrieve(req->datatype), req->count);
                }
                interposer_buffer_free(req->tmpbuf);
            }


            g_outstanding_requests.erase(req);
            break;
        }
    }
}

//**********************************************************


//TODO we should also implement testany, testsome, waitany, waitsome to be complete, but who uses those :D
int MPI_Wait(MPI_Request *request, MPI_Status *status) {
    int ret;

    if (*request != MPI_REQUEST_NULL) {
        ret = PMPI_Wait(request, status);
        interposer_request_free(request);
    }
    else {
        ret = PMPI_Wait(request, status);
    }

    return ret;
}

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
            interposer_request_free(&array_of_requests[i]);
        }

        /*
             //check if it was a Request, if yes, unpack TODO this is expensive :-(
            std::map<MPI_Request, struct Request>::iterator it;
            it = g_my_recv_requests.find(oldrequests[i]);
            if (it != g_my_recv_requests.end()) {
                DDT_Unpack(it->second.tmpbuf, it->second.usrbuf, datatype_retrieve(it->second.datatype), it->second.count);
                g_my_recv_requests.erase(it);
            }
            interposer_buffer_free(g_my_buffers[oldrequests[i]]);
            g_my_buffers.erase(oldrequests[i]);
        }
        */
    }

    free(oldrequests);

    return ret;
}


// TODO: MPI_Test and MPI_TestAll look buggy to me... How can we free buffers if
//       we don't check whether the test passed?
int MPI_Test(MPI_Request *request, int *flag, MPI_Status *status) {
    printf("Whoaa!");

    int ret;

    if (*request != MPI_REQUEST_NULL) {
        ret = PMPI_Test(request, flag, status);
        interposer_request_free(request);
        
        /*    
        The following line must surely be a bug?
//        if (*request == MPI_REQUEST_NULL) {
             //check if it was a Request, if yes, unpack TODO this is expensive :-(
            std::map<MPI_Request, struct Request>::iterator it;
            it = g_my_recv_requests.find(oldrequest);
            if (it != g_my_recv_requests.end()) {
                DDT_Unpack(it->second.tmpbuf, it->second.usrbuf, datatype_retrieve(it->second.datatype), it->second.count);
                g_my_recv_requests.erase(it);
            }
           
            interposer_buffer_free(g_my_buffers[oldrequest]);
            g_my_buffers.erase(oldrequest);
        }
        */
    }
    else {
        ret = PMPI_Test(request, flag, status);
    }

    return ret;

}

int MPI_Testall(int count, MPI_Request *array_of_requests, int *flag, MPI_Status *array_of_statuses) {
    printf("Whoaa!");

    // we need to copy the old requests here :-( (they will become MPI_REQUEST_NULL)
    MPI_Request* oldrequests = (MPI_Request*) malloc(count * sizeof(MPI_Request));
    for (int i=0; i<count; i++) {
        oldrequests[i] = array_of_requests[i];
    }

    int ret = PMPI_Testall(count, array_of_requests, flag, array_of_statuses);

    // now go over them and unpack recv requests and free buffers
    for (int i=0; i<count; i++) {
        if ((oldrequests[i] != MPI_REQUEST_NULL) && (array_of_requests[i] == MPI_REQUEST_NULL)) {

            interposer_request_free(&array_of_requests[i]);

            /*
             //check if it was a Request, if yes, unpack TODO this is expensive :-(
            std::map<MPI_Request, struct Request>::iterator it;
            it = g_my_recv_requests.find(oldrequests[i]);
            if (it != g_my_recv_requests.end()) {
                DDT_Unpack( it->second.tmpbuf, it->second.usrbuf, datatype_retrieve(it->second.datatype), it->second.count);
                g_my_recv_requests.erase(it);
            }
            interposer_buffer_free(g_my_buffers[oldrequests[i]]);
            g_my_buffers.erase(oldrequests[i]);
            */

        }
    }

    free(oldrequests);

    return ret;

}




// Fortran Bindings
extern "C" {
    void interposer_init_() {
        interposer_init();
    }

    void interposer_finalize_() {
        interposer_finalize();
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

    void* interposer_buffer_alloc_(int count, MPI_Datatype datatype, int* buf_size) {
        return interposer_buffer_alloc(count, datatype, buf_size);
    }

    void interposer_buffer_free_(void* tmpbuf) {
    printf("There\n");
        interposer_buffer_free(tmpbuf);
    }

    void interposer_request_register_(void *tmpbuf, void *usrbuf, int count, MPI_Datatype datatype, void* buf, MPI_Request *request) {
        interposer_request_register(tmpbuf, usrbuf, count, datatype, request);
    }

    void* interposer_pack_(void *data, int count, MPI_Datatype datatype, int *buf_size) {
        return interposer_pack(data, count, datatype, buf_size);
    }

    void interposer_unpack_(void *data, int count, MPI_Datatype datatype, void* buf) {
        interposer_unpack(data, count, datatype, buf);
    }
}
