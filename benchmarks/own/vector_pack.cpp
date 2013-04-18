// Copyright 2013 Timo Schneider and Fredrik Berg Kjolstad
//
// This file is part of the libpack packing library.
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT in the top level directory for details.

#include <string>
#include <mpi.h>
#include <vector>
#include <algorithm>

#include <ddt_jit.hpp>
#include "../../tests/test.hpp"
#include "../../copy_benchmark/hrtimer/hrtimer.h"


unsigned long long g_timerfreq;

void init_in_and_out_buffer(size_t buffer_size, char** inbuf, char** outbuf) {

	posix_memalign(reinterpret_cast<void**>(inbuf), 16, buffer_size);
  posix_memalign(reinterpret_cast<void**>(outbuf), 16, buffer_size);

	assert(inbuf != NULL);
	assert(outbuf != NULL);

    for (size_t i=0; i<buffer_size; i++) {
        ((char*) *inbuf)[i] = i+1;
        ((char*) *outbuf)[i] = 0;
    }
	
}

void free_in_and_out_buffer(size_t buffer_size, char** inbuf, char** outbuf) {

	free(*inbuf);
	free(*outbuf);
	
}

void benchmark_vector(int blklen, int stride, int inner_cnt, int outer_cnt, int inner_runs, int outer_runs) {

    char* farc_inbuf;
    char* farc_outbuf;

    farc::Datatype* t1 = new farc::PrimitiveDatatype(farc::PrimitiveDatatype::DOUBLE);
    farc::Datatype* t2 = new farc::VectorDatatype(t1, inner_cnt, blklen, stride);
    farc::DDT_Commit(t2);
    int data_size   = t2->getSize();
    int buffer_size = t2->getExtent();

    init_in_and_out_buffer(buffer_size, &farc_inbuf, &farc_outbuf);

    farc::DDT_Pack(farc_inbuf, farc_outbuf, t2, outer_cnt);

    farc::DDT_Free(t1);
    farc::DDT_Free(t2);

    free_in_and_out_buffer(buffer_size, &farc_inbuf, &farc_outbuf);

}

int main(int argc, char** argv) {

   if (argc < 4) {
        fprintf(stderr, "%s [outer_cnt] [count] [blklen] [stride]\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    farc::DDT_Init();
 
    int inner_cnt = atoi(argv[1]);
    int blklen = atoi(argv[2]);
    int stride = atoi(argv[3]);

    benchmark_vector(blklen, stride, inner_cnt, 1, 1, 1);

    return 0;

}

