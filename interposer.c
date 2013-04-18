// Copyright 2013 Timo Schneider and Fredrik Berg Kjolstad
//
// This file is part of the libpack packing library.
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT in the top level directory for details.

/**
Link this into your application to use FARC DDTs (see the pmpi-tests Makefile for how to do that)
*/

#include <mpi.h>
#include <stdlib.h>

#include "interposer_common.h"

int MPI_Init(int *argc, char ***argv) {
    interposer_init();
    return PMPI_Init(argc, argv);
}

int MPI_Finalize(void) {
    interposer_finalize();
    return PMPI_Finalize();
}

int MPI_Type_hvector(int count, int blocklength, MPI_Aint stride, MPI_Datatype oldtype, MPI_Datatype *newtype) {
//    PMPI_Type_hvector(count, blocklength, stride, oldtype, newtype);
    interposer_hvector(count, blocklength, stride, oldtype, newtype);
    return MPI_SUCCESS;
}

int MPI_Type_vector(int count, int blocklength, int stride, MPI_Datatype oldtype, MPI_Datatype *newtype) {
//    PMPI_Type_vector(count, blocklength, stride, oldtype, newtype);
    interposer_vector(count, blocklength, stride, oldtype, newtype);
    return MPI_SUCCESS;
}

int MPI_Type_create_struct(int count, int array_of_blocklengths[], MPI_Aint array_of_displacements[], MPI_Datatype array_of_types[], MPI_Datatype *newtype) {
    interposer_create_struct(count, array_of_blocklengths, array_of_displacements, array_of_types, newtype);
    return MPI_SUCCESS;
}

int MPI_Type_struct(int count, int *array_of_blocklengths, MPI_Aint *array_of_displacements, MPI_Datatype *array_of_types, MPI_Datatype *newtype) {
    interposer_struct(count, array_of_blocklengths, array_of_displacements, array_of_types, newtype);
    return MPI_SUCCESS;
}

int MPI_Type_create_hvector(int count, int blocklength, MPI_Aint stride, MPI_Datatype oldtype, MPI_Datatype *newtype) {
    interposer_create_hvector(count, blocklength, stride, oldtype, newtype);
    return MPI_SUCCESS;
}

int MPI_Type_create_hindexed(int count, int array_of_blocklengths[], MPI_Aint array_of_displacements[], MPI_Datatype oldtype, MPI_Datatype *newtype) {
    interposer_create_hindexed(count, array_of_blocklengths, array_of_displacements, oldtype, newtype);
    return MPI_SUCCESS;
}

int MPI_Type_hindexed(int count, int array_of_blocklengths[], MPI_Aint array_of_displacements[], MPI_Datatype oldtype, MPI_Datatype *newtype) {
    interposer_hindexed(count, array_of_blocklengths, array_of_displacements, oldtype, newtype);
    return MPI_SUCCESS;
}

int MPI_Type_create_indexed_block(int count, int blocklength, int array_of_displacements[], MPI_Datatype oldtype, MPI_Datatype *newtype) {
    interposer_indexed_block(count, blocklength, array_of_displacements, oldtype, newtype);
    return MPI_SUCCESS;    
}

int MPI_Type_contiguous(int count, MPI_Datatype oldtype, MPI_Datatype *newtype) {
    interposer_contiguous(count, oldtype, newtype);
    return MPI_SUCCESS;
}

int MPI_Type_commit(MPI_Datatype *datatype) {
    interposer_commit(datatype);
    return MPI_SUCCESS;
}

int MPI_Type_free(MPI_Datatype *datatype) {
    interposer_free(datatype);
    return MPI_SUCCESS;
}

int MPI_Type_size(MPI_Datatype datatype, int *size) {

    *size = interposer_type_size(datatype);

    return MPI_SUCCESS;

}

int MPI_Type_extent(MPI_Datatype datatype, int *extent) {

    *extent = interposer_type_extent(datatype);

    return MPI_SUCCESS;

}

int MPI_Send(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm) {
    // TODO check other primitive types!
    if ((datatype != MPI_DOUBLE) && (datatype != MPI_FLOAT) && (datatype != MPI_INT) && (datatype != MPI_BYTE) && (datatype != MPI_CHAR) && (datatype != MPI_PACKED)) {
        int outsize;
        void *outbuf = interposer_pack(buf, count, datatype, &outsize);

        PMPI_Send(outbuf, outsize, MPI_BYTE, dest, tag, comm);

        interposer_buffer_free(outbuf);
        return MPI_SUCCESS;
    }
    else {
        PMPI_Send(buf, count, datatype, dest, tag, comm);
    }

    return MPI_SUCCESS;
}

int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status) {
    // TODO check other primitive types!
    if ((datatype != MPI_DOUBLE) && (datatype != MPI_FLOAT) && (datatype != MPI_INT) && (datatype != MPI_BYTE) && (datatype != MPI_CHAR) && (datatype != MPI_PACKED)) {
        int insize;
        void* inbuf = interposer_buffer_alloc(count, datatype, &insize);

        PMPI_Recv(inbuf, insize, MPI_BYTE, source, tag, comm, status);

        interposer_unpack(buf, count, datatype, inbuf);
        interposer_buffer_free(inbuf);

        return MPI_SUCCESS;
    }
    else {
        PMPI_Recv(buf, count, datatype, source, tag, comm, status);
    }

    return MPI_SUCCESS;
}

int MPI_Isend(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request) {
    // TODO check other primitive types!
    if ((datatype != MPI_DOUBLE) && (datatype != MPI_FLOAT) && (datatype != MPI_INT) && (datatype != MPI_BYTE) && (datatype != MPI_CHAR) && (datatype != MPI_PACKED)) {
        int outsize;
        void *outbuf = interposer_pack(buf, count, datatype, &outsize);
    
        PMPI_Isend(outbuf, outsize, MPI_BYTE, dest, tag, comm, request);

        interposer_request_register(outbuf, NULL, 0, 0, request);
        return MPI_SUCCESS;
    }
    else {
        PMPI_Isend(buf, count, datatype, dest, tag, comm, request);
        interposer_request_register(NULL, NULL, 0, 0, request);
    }

    return MPI_SUCCESS;
}

int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request) {
    // TODO check other primitive types!
    if ((datatype != MPI_DOUBLE) && (datatype != MPI_FLOAT) && (datatype != MPI_INT) && (datatype != MPI_BYTE) && (datatype != MPI_CHAR) && (datatype != MPI_PACKED)) {
        int insize;
        void* inbuf = interposer_buffer_alloc(count, datatype, &insize);

        PMPI_Irecv(inbuf, insize, MPI_BYTE, source, tag, comm, request);

        interposer_request_register(inbuf, buf, count, datatype, request);
        return MPI_SUCCESS;
    }
    else {
        PMPI_Irecv(buf, count, datatype, source, tag, comm, request);
        interposer_request_register(NULL, NULL, 0, 0, request);
    }

    return MPI_SUCCESS;
  
}

int MPI_Pack(void *inbuf, int incount, MPI_Datatype datatype, void *outbuf, int outsize, int *position, MPI_Comm comm) {

    //TODO we ignore position and outsize
    interposer_pack_providedbuf(inbuf, incount, datatype, outbuf);

    return MPI_SUCCESS;

}

int MPI_Unpack(void *inbuf, int insize, int *position, void *outbuf, int outcount, MPI_Datatype datatype, MPI_Comm comm) {

    //TODO we ignore position
    interposer_unpack(outbuf, outcount, datatype, inbuf);

    return MPI_SUCCESS;

}


