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

    char* mpi_inbuf;
    char* mpi_outbuf;
    char* farc_inbuf;
    char* farc_outbuf;
    double* cpp_inbuf;
    double* cpp_outbuf;
    int jend, j;

    HRT_TIMESTAMP_T start, stop;
    std::vector<uint64_t> mpi_type_create (inner_runs*outer_runs, 0);
	std::vector<uint64_t> farc_type_create (inner_runs*outer_runs, 0);
	std::vector<uint64_t> mpi_pack (inner_runs*outer_runs, 0);
	std::vector<uint64_t> farc_pack (inner_runs*outer_runs, 0);
	std::vector<uint64_t> cpp_pack (inner_runs*outer_runs, 0);

    farc::Datatype* tmp1 = new farc::PrimitiveDatatype(farc::PrimitiveDatatype::DOUBLE);
    farc::Datatype* tmp2 = new farc::VectorDatatype(tmp1, inner_cnt, blklen, stride);
    int data_size   = tmp2->getSize();
    int buffer_size = tmp2->getExtent();
    farc::DDT_Free(tmp1);
    farc::DDT_Free(tmp2);

    init_in_and_out_buffer(buffer_size, &mpi_inbuf, &mpi_outbuf);
    init_in_and_out_buffer(buffer_size, &farc_inbuf, &farc_outbuf);
    init_in_and_out_buffer(buffer_size, reinterpret_cast<char**>(&cpp_inbuf),  reinterpret_cast<char**>(&cpp_outbuf));

    for (int o=0; o<outer_runs; o++) {
        HRT_GET_TIMESTAMP(start);
        farc::Datatype* t1 = new farc::PrimitiveDatatype(farc::PrimitiveDatatype::DOUBLE);
        farc::Datatype* t2 = new farc::VectorDatatype(t1, inner_cnt, blklen, stride);
        farc::DDT_Commit(t2);
        HRT_GET_TIMESTAMP(stop);
        HRT_GET_ELAPSED_TICKS(start, stop, &farc_type_create[o*inner_runs]);

        HRT_GET_TIMESTAMP(start);
        MPI_Datatype newtype_mpi;
        MPI_Type_vector(inner_cnt, blklen, stride, MPI_DOUBLE, &newtype_mpi);
        MPI_Type_commit(&newtype_mpi);
        HRT_GET_TIMESTAMP(stop);
        HRT_GET_ELAPSED_TICKS(start, stop, &mpi_type_create[o*inner_runs]);

        for (int i=0; i<inner_runs; i++) {
            int position = 0;

			//warmup
            MPI_Pack(mpi_inbuf, outer_cnt, newtype_mpi, mpi_outbuf, buffer_size*sizeof(int), &position, MPI_COMM_WORLD);

            HRT_GET_TIMESTAMP(start);
            MPI_Pack(mpi_inbuf, outer_cnt, newtype_mpi, mpi_outbuf, buffer_size*sizeof(int), &position, MPI_COMM_WORLD);
            HRT_GET_TIMESTAMP(stop);
            HRT_GET_ELAPSED_TICKS(start, stop, &mpi_pack[i+o*inner_runs]);

			// warmup
            for (int ii=0; ii < outer_cnt; ii++) {
                for(int jj=0; jj < inner_cnt; jj++) {
                    for(int kk=0; kk < blklen; kk++) {
                        cpp_outbuf[ii*inner_cnt*blklen + jj*blklen + kk] = cpp_inbuf[ii*((inner_cnt-1)*stride + blklen) + jj*stride + kk];
                    }
                }
            }


            HRT_GET_TIMESTAMP(start);
            for (int ii=0; ii < outer_cnt; ii++) {
                for(int jj=0; jj < inner_cnt; jj++) {
                    for(int kk=0; kk < blklen; kk++) {
                        cpp_outbuf[ii*inner_cnt*blklen + jj*blklen + kk] = cpp_inbuf[ii*((inner_cnt-1)*stride + blklen) + jj*stride + kk];
                    }
                }
            }
            HRT_GET_TIMESTAMP(stop);
            HRT_GET_ELAPSED_TICKS(start, stop, &cpp_pack[i+o*inner_runs]);

			// warmup
            farc::DDT_Pack(farc_inbuf, farc_outbuf, t2, outer_cnt);

            HRT_GET_TIMESTAMP(start);
            farc::DDT_Pack(farc_inbuf, farc_outbuf, t2, outer_cnt);
            HRT_GET_TIMESTAMP(stop);
            HRT_GET_ELAPSED_TICKS(start, stop, &farc_pack[i+o*inner_runs]);

			mpi_type_create[i+o*inner_runs] = mpi_type_create[o*inner_runs];
			farc_type_create[i+o*inner_runs] = farc_type_create[o*inner_runs];

        } 
        MPI_Type_free(&newtype_mpi);
        farc::DDT_Free(t1);
        farc::DDT_Free(t2);

	  }

	  std::sort(mpi_type_create.begin(), mpi_type_create.end());
	  std::sort(farc_type_create.begin(), farc_type_create.end());
	  std::sort(cpp_pack.begin(), cpp_pack.end());
	  std::sort(mpi_pack.begin(), mpi_pack.end());
	  std::sort(farc_pack.begin(), farc_pack.end());


	  int r = inner_runs * outer_runs;
      static int firstline=1;
      if (firstline) printf("%10s %10s %10s %10s %10s %10s %10s %10s %10s %10s \n", "size", "mpi_create", "farc_create", "cpp_pack", "mpi_pack", "farc_pack", "blklen", "stride", "count", "pack_count");
      firstline=0;
      printf("%10i %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10i %10i %10i %10i\n", data_size, HRT_GET_USEC(mpi_type_create[r/2]), 
	  																					      HRT_GET_USEC(farc_type_create[r/2]), 
																							  HRT_GET_USEC(cpp_pack[r/2]), 
																							  HRT_GET_USEC(mpi_pack[r/2]), 
																							  HRT_GET_USEC(farc_pack[r/2]), blklen, stride, inner_cnt, outer_cnt);



	free_in_and_out_buffer(buffer_size, &mpi_inbuf, &mpi_outbuf);
	free_in_and_out_buffer(buffer_size, &farc_inbuf, &farc_outbuf);
	free_in_and_out_buffer(buffer_size, reinterpret_cast<char**>(&cpp_inbuf), reinterpret_cast<char**>(&cpp_outbuf));

}

int main(int argc, char** argv) {

   if (argc < 14) {
        fprintf(stderr, "%s [num-runs] [blklen_start] [blklen_end] [blklen_inc] [stride_start] [stride_end] [stride_inc] [inner_cnt_start] [inner_cnt_end] [inner_cnt_inc] [outer_cnt_start] [outer_cnt_end] [outer_cnt_inc]\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    MPI_Init(&argc, &argv);
    HRT_INIT(1, g_timerfreq);
    farc::DDT_Init();
 
    int num_runs =  atoi(argv[1]);
    srand(time(NULL));

    int blklen_start = atoi(argv[2]);
    int blklen_end   = atoi(argv[3]);
    int blklen_inc   = atoi(argv[4]);

    int stride_start = atoi(argv[5]);
    int stride_end   = atoi(argv[6]);
    int stride_inc   = atoi(argv[7]);

    int inner_cnt_start = atoi(argv[8]);
    int inner_cnt_end   = atoi(argv[9]);
    int inner_cnt_inc   = atoi(argv[10]);

    int outer_cnt_start = atoi(argv[11]);
    int outer_cnt_end = atoi(argv[12]);
    int outer_cnt_inc = atoi(argv[13]);

    for (int blklen=blklen_start; blklen<=blklen_end; blklen += blklen_inc) {
        for (int stride=stride_start; stride<=stride_end; stride += stride_inc) {
            for (int inner_cnt=inner_cnt_start; inner_cnt<=inner_cnt_end; inner_cnt += inner_cnt_inc) {
                for (int outer_cnt=outer_cnt_start; outer_cnt<=outer_cnt_end; outer_cnt += outer_cnt_inc) {
                    benchmark_vector(blklen, stride, inner_cnt, outer_cnt, 5, num_runs);
                }
            }
        }
    }

    MPI_Finalize();
    return 0;

}

