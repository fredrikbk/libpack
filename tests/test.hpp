#ifndef TEST_HPP
#define TEST_HPP

#include <assert.h>
#include <cstdio>
#include <cstdlib>
#include <string>

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

int compare_buffers(size_t size, char** mpi_inbuf, char** farc_inbuf, char** mpi_outbuf, char** farc_outbuf) {

    for (size_t i=0; i<size; i++) {
        if ((*mpi_inbuf)[i] != (*farc_inbuf)[i]) return -1;
        if ((*mpi_outbuf)[i] != (*farc_outbuf)[i]) return -1;
    }   

    return 0; 

}

int inspect_buffers(size_t size, char** mpi_inbuf, char** farc_inbuf, char** mpi_outbuf, char** farc_outbuf) {

    for (int i=0; i<size; i++) {
        printf("mpi_inbuf[%i] = %i farc_inbuf[%i] = %i mpi_outbuf[%i] = %i farc_outbuf[%i] = %i\n", i, (*mpi_inbuf)[i], i, (*farc_inbuf)[i], i, (*mpi_outbuf)[i], i, (*farc_outbuf)[i]);
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

void free_buffers(char** mpi_inbuf, char** farc_inbuf, char** mpi_outbuf, char** farc_outbuf) {

    free(*mpi_inbuf);
    free(*farc_inbuf);
    free(*mpi_outbuf);
    free(*farc_outbuf);

}

#endif

