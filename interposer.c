/**
Link this into your application to use FARC DDTs (see the pmpi-tests Makefile for how to do that)
*/

#include <mpi.h>
#include "interposer_common.h"

int MPI_Init(int *argc, char ***argv) {
    interposer_init();
    return PMPI_Init(argc, argv);
}

int MPI_Type_hvector(int count, int blocklength, MPI_Aint stride, MPI_Datatype oldtype, MPI_Datatype *newtype) {
    PMPI_Type_hvector(count, blocklength, stride, oldtype, newtype);
    interposer_hvector(count, blocklength, stride, oldtype, newtype);
    return MPI_SUCCESS;
}

int MPI_Type_vector(int count, int blocklength, int stride, MPI_Datatype oldtype, MPI_Datatype *newtype) {
    PMPI_Type_vector(count, blocklength, stride, oldtype, newtype);
    interposer_vector(count, blocklength, stride, oldtype, newtype);
    return MPI_SUCCESS;
}

int MPI_Type_create_struct(int count, int array_of_blocklengths[], MPI_Aint array_of_displacements[], MPI_Datatype array_of_types[], MPI_Datatype *newtype) {
    PMPI_Type_create_struct(count, array_of_blocklengths, array_of_displacements, array_of_types, newtype);
    interposer_create_struct(count, array_of_blocklengths, array_of_displacements, array_of_types, newtype);
    return MPI_SUCCESS;
}

int MPI_Type_struct(int count, int *array_of_blocklengths, MPI_Aint *array_of_displacements, MPI_Datatype *array_of_types, MPI_Datatype *newtype) {
    PMPI_Type_struct(count, array_of_blocklengths, array_of_displacements, array_of_types, newtype);
    interposer_struct(count, array_of_blocklengths, array_of_displacements, array_of_types, newtype);
    return MPI_SUCCESS;
}

int MPI_Type_create_hvector(int count, int blocklength, MPI_Aint stride, MPI_Datatype oldtype, MPI_Datatype *newtype) {
    PMPI_Type_create_hvector(count, blocklength, stride, oldtype, newtype);
    interposer_create_hvector(count, blocklength, stride, oldtype, newtype);
    return MPI_SUCCESS;
}

int MPI_Type_create_hindexed(int count, int array_of_blocklengths[], MPI_Aint array_of_displacements[], MPI_Datatype oldtype, MPI_Datatype *newtype) {
    PMPI_Type_create_hindexed(count, array_of_blocklengths, array_of_displacements, oldtype, newtype);
    interposer_create_hindexed(count, array_of_blocklengths, array_of_displacements, oldtype, newtype);
    return MPI_SUCCESS;
}

int MPI_Type_hindexed(int count, int array_of_blocklengths[], MPI_Aint array_of_displacements[], MPI_Datatype oldtype, MPI_Datatype *newtype) {
    PMPI_Type_hindexed(count, array_of_blocklengths, array_of_displacements, oldtype, newtype);
    interposer_hindexed(count, array_of_blocklengths, array_of_displacements, oldtype, newtype);
    return MPI_SUCCESS;
}

int MPI_Type_contiguous(int count, MPI_Datatype oldtype, MPI_Datatype *newtype) {
    PMPI_Type_contiguous(count, oldtype, newtype);
    interposer_contiguous(count, oldtype, newtype);
    return MPI_SUCCESS;
}

int MPI_Type_commit(MPI_Datatype *datatype) {
    PMPI_Type_commit(datatype);
    interposer_commit(datatype);
    return MPI_SUCCESS;
}

int MPI_Type_free(MPI_Datatype *datatype) {
    interposer_free(datatype);
    return PMPI_Type_free(datatype);
}

int MPI_Send(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm) {
    // TODO check other primitive types!
    if (datatype != MPI_DOUBLE) {
	int outsize;
	char *outbuf = interposer_pack(buf, count, datatype, &outsize);

        PMPI_Send(outbuf, outsize, MPI_BYTE, dest, tag, comm);

        interposer_buffer_free(outbuf);
        return MPI_SUCCESS;
    }
    else {
        PMPI_Send(buf, count, MPI_DOUBLE, dest, tag, comm);
    }

    return MPI_SUCCESS;
}

int MPI_Isend(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request) {
    // TODO check other primitive types!
    if (datatype != MPI_DOUBLE) {
	int outsize;
	char *outbuf = interposer_pack(buf, count, datatype, &outsize);
    
        PMPI_Isend(outbuf, outsize, MPI_BYTE, dest, tag, comm, request);

	interposer_buffer_register(request, outbuf);
        return MPI_SUCCESS;
    }
    else {
        PMPI_Isend(buf, count, MPI_DOUBLE, dest, tag, comm, request);
	interposer_buffer_register(request, NULL);
    }

    return MPI_SUCCESS;
}

int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status) {
    // TODO check other primitive types!
    if (datatype != MPI_DOUBLE) {
	int insize;
	char* inbuf = interposer_buffer_alloc(count, datatype, &insize);

        PMPI_Recv(inbuf, insize, MPI_BYTE, source, tag, comm, status);

	interposer_unpack(buf, count, datatype, inbuf);
	interposer_buffer_free(inbuf);

        return MPI_SUCCESS;
    }
    else {
        PMPI_Recv(buf, count, MPI_DOUBLE, source, tag, comm, status);
    }

    return MPI_SUCCESS;
}

int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request) {

    // TODO check other primitive types!

    if (datatype != MPI_DOUBLE) {
	int insize;
	char* inbuf = interposer_buffer_alloc(count, datatype, &insize);

        PMPI_Irecv(inbuf, insize, MPI_BYTE, source, tag, comm, request);

	interposer_buffer_register(request, inbuf);
	interposer_recvop_register(buf, count, datatype, inbuf, request);

        return MPI_SUCCESS;
    }
    else {
        PMPI_Irecv(buf, count, MPI_DOUBLE, source, tag, comm, request);
	interposer_buffer_register(request, NULL);
    }

    return MPI_SUCCESS;
  
}

