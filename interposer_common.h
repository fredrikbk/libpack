#ifndef INTERPOSER_COMMON_H 
#define INTERPOSER_COMMON_H 

#include <mpi.h>

void interposer_init();
void interposer_finalize();
void interposer_hvector(int count, int blocklength, MPI_Aint stride, MPI_Datatype oldtype, MPI_Datatype *newtype);
void interposer_vector(int count, int blocklength, int stride, MPI_Datatype oldtype, MPI_Datatype *newtype);
void interposer_create_struct(int count, int array_of_blocklengths[], MPI_Aint array_of_displacements[], MPI_Datatype array_of_types[], MPI_Datatype *newtype);
void interposer_struct(int count, int *array_of_blocklengths, MPI_Aint *array_of_displacements, MPI_Datatype *array_of_types, MPI_Datatype *newtype);
void interposer_create_hvector(int count, int blocklength, MPI_Aint stride, MPI_Datatype oldtype, MPI_Datatype *newtype);
void interposer_create_hindexed(int count, int array_of_blocklengths[], MPI_Aint array_of_displacements[], MPI_Datatype oldtype, MPI_Datatype *newtype);
void interposer_hindexed(int count, int array_of_blocklengths[], MPI_Aint array_of_displacements[], MPI_Datatype oldtype, MPI_Datatype *newtype);
void interposer_indexed_block(int count, int blocklength, int array_of_displacements[], MPI_Datatype oldtype, MPI_Datatype *newtype);
void interposer_contiguous(int count, MPI_Datatype oldtype, MPI_Datatype *newtype);
void interposer_commit(MPI_Datatype *datatype);
void interposer_free(MPI_Datatype *datatype);
int interposer_type_size(MPI_Datatype datatype);
int interposer_type_extent(MPI_Datatype datatype);
void* interposer_buffer_alloc(int count, MPI_Datatype datatype, int* buf_size);
void interposer_buffer_free(void* buf);
void interposer_request_register(void *tmpbuf, void *usrbuf, int count, MPI_Datatype datatype, MPI_Request *request);
void* interposer_pack(void *data, int count, MPI_Datatype datatype, int *buf_size);
void interposer_pack_providedbuf(void* inbuf, int incount, MPI_Datatype datatype, void *outbuf);
void interposer_unpack(void *data, int count, MPI_Datatype datatype, void* buf);

#endif
