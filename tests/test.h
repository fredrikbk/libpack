// Copyright 2013 Timo Schneider and Fredrik Berg Kjolstad
//
// This file is part of the libpack packing library.
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT in the top level directory for details.

#ifndef TEST_HPP
#define TEST_HPP

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

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

void init_buffers_aligned(size_t size, char** mpi_inbuf, char** farc_inbuf, char** mpi_outbuf, char** farc_outbuf) {

    srand(time(NULL));
    
    int ret;
    ret = posix_memalign((void**) mpi_inbuf, 16, size);
    assert(ret == 0);
    ret = posix_memalign((void**) mpi_outbuf, 16, size);
    assert(ret == 0);
    ret = posix_memalign((void**) farc_inbuf, 16, size);
    assert(ret == 0);
    ret = posix_memalign((void**) farc_outbuf, 16, size);
    assert(ret == 0);

    size_t i;
    for (i=0; i<size; i++) {
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

    size_t i;
    for (i=0; i<size; i++) {
        printf("mpi_inbuf[%i] = %i farc_inbuf[%i] = %i mpi_outbuf[%i] = %i farc_outbuf[%i] = %i (@%p)\n", 
                i, (*mpi_inbuf)[i], i, (*farc_inbuf)[i], i, (*mpi_outbuf)[i], i, (*farc_outbuf)[i], &((*farc_outbuf)[i]));
    }   

    return 0; 

}

void test_start(char* str) {
    printf("%100s", str);
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

#endif

