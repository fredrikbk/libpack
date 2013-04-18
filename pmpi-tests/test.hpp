// Copyright 2013 Timo Schneider and Fredrik Berg Kjolstad
//
// This file is part of the libpack packing library.
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT in the top level directory for details.

#ifndef TEST_HPP
#define TEST_HPP


#include <mpi.h>
#include <assert.h>

#ifdef __cplusplus
#include <cstdio>
#include <cstdlib>
#include <string>
#else
#include <stdio.h>
#include <stdlib.h>
#endif

void init_buffers(size_t size, char** mpi_inbuf, char** farc_inbuf, char** mpi_outbuf, char** farc_outbuf) {

    srand(time(NULL));
    
    *mpi_inbuf = (char*) malloc(size);
    *mpi_outbuf = (char*) malloc(size);
    *farc_inbuf = (char*) malloc(size);
    *farc_outbuf = (char*) malloc(size);
    assert(*mpi_inbuf != NULL);
    assert(*farc_inbuf != NULL);

    size_t i;
    for (i=0; i<size; i++) {
        (*mpi_inbuf)[i] = i+1;
        (*farc_inbuf)[i] = i+1;
        (*mpi_outbuf)[i] = 0;
        (*farc_outbuf)[i] = 0;
    }

}

int compare_buffers(size_t size, char** mpi_inbuf, char** farc_inbuf, char** mpi_outbuf, char** farc_outbuf) {

    size_t i;
    for (i=0; i<size; i++) {
        if ((*mpi_inbuf)[i] != (*farc_inbuf)[i]) return -1;
        if ((*mpi_outbuf)[i] != (*farc_outbuf)[i]) return -1;
    }   

    return 0; 

}

int inspect_buffers(size_t size, char** mpi_inbuf, char** farc_inbuf, char** mpi_outbuf, char** farc_outbuf) {

    printf("\n");
    int i;
    for (i=0; i<size; i++) {
        printf("mpi_inbuf[%i] = %i farc_inbuf[%i] = %i mpi_outbuf[%i] = %i farc_outbuf[%i] = %i\n", i, (*mpi_inbuf)[i], i, (*farc_inbuf)[i], i, (*mpi_outbuf)[i], i, (*farc_outbuf)[i]);
    }   

    return 0; 

}

#ifdef __cplusplus
void test_start(std::string str) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) printf("%100s", str.c_str());

}
#else
void test_start(char* str) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) printf("%100s", str);
}
#endif

void test_result(int res) {
    int rank, res_global;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Reduce(&res, &res_global, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    if ((rank == 0) && (res_global == 0)) printf(" ... ok\n");
    else if (rank == 0) printf(" ... not ok\n");
}

void free_buffers(char** mpi_inbuf, char** farc_inbuf, char** mpi_outbuf, char** farc_outbuf) {

    free(*mpi_inbuf);
    free(*farc_inbuf);
    free(*mpi_outbuf);
    free(*farc_outbuf);

}

#endif

