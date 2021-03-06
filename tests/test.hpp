// Copyright 2013 Timo Schneider and Fredrik Berg Kjolstad
//
// This file is part of the libpack packing library.
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT in the top level directory for details.

#ifndef TEST_HPP
#define TEST_HPP

#include <assert.h>
#include <cstdio>
#include <cstdlib>
#include <string>
#include "../ddt_jit.hpp"

void init_buffers(size_t size, char** mpi_inbuf, char** farc_inbuf, char** mpi_outbuf, char** farc_outbuf) {

    srand(time(NULL));
    
    *mpi_inbuf = (char*) malloc(size);
    *mpi_outbuf = (char*) malloc(size);
    *farc_inbuf = (char*) malloc(size);
    *farc_outbuf = (char*) malloc(size);
    assert(*mpi_inbuf != NULL);
    assert(*farc_inbuf != NULL);

    for (size_t i=0; i<size; i++) {
        (*mpi_inbuf)[i] = i+1;
        (*farc_inbuf)[i] = i+1;
        (*mpi_outbuf)[i] = 0;
        (*farc_outbuf)[i] = 0;
    }

}

void init_buffers_aligned(size_t size, char** mpi_inbuf, char** farc_inbuf, char** mpi_outbuf, char** farc_outbuf) {

    srand(time(NULL));
    
    posix_memalign(reinterpret_cast<void**>(mpi_inbuf), 16, size);
    posix_memalign(reinterpret_cast<void**>(mpi_outbuf), 16, size);
    posix_memalign(reinterpret_cast<void**>(farc_inbuf), 16, size);
    posix_memalign(reinterpret_cast<void**>(farc_outbuf), 16, size);

    for (size_t i=0; i<size; i++) {
        (*mpi_inbuf)[i] = i+1;
        (*farc_inbuf)[i] = i+1;
        (*mpi_outbuf)[i] = 0;
        (*farc_outbuf)[i] = 0;
    }

}

void init_buffers_unaligned(size_t size, char** mpi_inbuf, char** farc_inbuf, char** mpi_outbuf, char** farc_outbuf) {

    srand(time(NULL));
    
    *mpi_inbuf = (char*) malloc(size);
    *mpi_outbuf = (char*) malloc(size);
    *farc_inbuf = (char*) malloc(size);
    *farc_outbuf = (char*) malloc(size);

    *mpi_inbuf += 3;
    *mpi_outbuf += 3;
    *farc_inbuf += 3;
    *farc_outbuf += 3;

    assert(*mpi_inbuf != NULL);
    assert(*farc_inbuf != NULL);

    for (size_t i=0; i<size; i++) {
        (*mpi_inbuf)[i] = i+1;
        (*farc_inbuf)[i] = i+1;
        (*mpi_outbuf)[i] = 0;
        (*farc_outbuf)[i] = 0;
    }

}

int compare_buffers(size_t size, char** mpi_inbuf, char** farc_inbuf, char** mpi_outbuf, char** farc_outbuf) {

    for (size_t i=0; i<size; i++) {
        if ((*mpi_inbuf)[i] != (*farc_inbuf)[i]) return -1;
        if ((*mpi_outbuf)[i] != (*farc_outbuf)[i]) return -1;
    }   

    return 0; 

}

int inspect_buffers(size_t size, char** mpi_inbuf, char** farc_inbuf, char** mpi_outbuf, char** farc_outbuf) {

    printf("\n");
    for (unsigned int i=0; i<size; i++) {
        printf("mpi_inbuf[%i] = %i farc_inbuf[%i] = %i mpi_outbuf[%i] = %i farc_outbuf[%i] = %i (@%p)\n", 
                i, (*mpi_inbuf)[i], i, (*farc_inbuf)[i], i, (*mpi_outbuf)[i], i, (*farc_outbuf)[i], &((*farc_outbuf)[i]));
    }   

    return 0; 

}

void test_start(std::string str) {
    printf("%100s", str.c_str());
}

void test_result(int res) {
    if (res == 0) printf(" ... ok\n");
    else printf(" ... not ok\n");
}

void free_buffers_unaligned(char** mpi_inbuf, char** farc_inbuf, char** mpi_outbuf, char** farc_outbuf) {

    free(*mpi_inbuf - 3);
    free(*farc_inbuf - 3);
    free(*mpi_outbuf - 3);
    free(*farc_outbuf - 3);

}

void free_buffers(char** mpi_inbuf, char** farc_inbuf, char** mpi_outbuf, char** farc_outbuf) {

    free(*mpi_inbuf);
    free(*farc_inbuf);
    free(*mpi_outbuf);
    free(*farc_outbuf);

}

int compare_ddt_info(MPI_Datatype mpitype, farc::Datatype* farctype) {

    MPI_Aint mpi_lb, mpi_extent;
    int farc_lb, farc_extent;

    MPI_Type_get_extent(mpitype, &mpi_lb, &mpi_extent);
    farc_lb = farctype->getLowerBound();
    farc_extent = farctype->getExtent();

    if ((farc_lb != mpi_lb) || (farc_extent != mpi_extent)) {
        return -1;
    }

    return 0;

}

void inspect_ddt_info(MPI_Datatype mpitype, farc::Datatype* farctype) {

    MPI_Aint mpi_lb, mpi_extent;
    int farc_lb, farc_extent;

    MPI_Type_get_extent(mpitype, &mpi_lb, &mpi_extent);
    farc_lb = farctype->getLowerBound();
    farc_extent = farctype->getExtent();

    printf("mpi extent: %i mpi lb: %i farc extent: %i farc lb: %i\n", 
                    mpi_extent, mpi_lb, farc_extent, farc_lb);

}


#endif

