// Copyright 2013 Timo Schneider and Fredrik Berg Kjolstad
//
// This file is part of the libpack packing library.
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT in the top level directory for details.

#include <string>
#include <mpi.h>

#include "../../ddt_jit.hpp"
#include "../../tests/test.hpp"
#include "../../copy_benchmark/hrtimer/hrtimer.h"

unsigned long long g_timerfreq;

void benchmark_vector(int blklen, int stride, int inner_cnt, int outer_cnt, int inner_runs, int outer_runs) {

    char* mpi_inbuf;
    char* mpi_outbuf;
    char* farc_inbuf;
    char* farc_outbuf;
    char* interposer_inbuf;
    char* interposer_outbuf;
    double* cpp_inbuf;
    double* cpp_outbuf;
    int jend, j;

    HRT_TIMESTAMP_T start, stop;
    uint64_t mpi_type_create, farc_type_create, interposer_type_create, mpi_send, farc_send, interposer_send, cpp_send;

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int data_size = sizeof(double)*inner_cnt * outer_cnt * blklen;
    int buffer_size = sizeof(double)*((inner_cnt-1)*stride+blklen) * outer_cnt;

    init_buffers(buffer_size, &mpi_inbuf, &farc_inbuf, &mpi_outbuf, &farc_outbuf);
    init_buffers(buffer_size, &interposer_inbuf, ((char**)&cpp_inbuf), &interposer_outbuf, ((char**)&cpp_outbuf));

    for (int o=0; o<outer_runs; o++) {
        HRT_GET_TIMESTAMP(start);
        farc::Datatype* t1 = new farc::PrimitiveDatatype(farc::PrimitiveDatatype::DOUBLE);
        farc::Datatype* t2 = new farc::VectorDatatype(t1, inner_cnt, blklen, stride);
        farc::DDT_Commit(t2);
        int t2_size = t2->getSize();
        HRT_GET_TIMESTAMP(stop);
        HRT_GET_ELAPSED_TICKS(start, stop, &farc_type_create);
	
        HRT_GET_TIMESTAMP(start);
        MPI_Datatype newtype;
        PMPI_Type_vector(inner_cnt, blklen, stride, MPI_DOUBLE, &newtype);
        PMPI_Type_commit(&newtype);
        HRT_GET_TIMESTAMP(stop);
        HRT_GET_ELAPSED_TICKS(start, stop, &mpi_type_create);

        HRT_GET_TIMESTAMP(start);
        MPI_Datatype newtype_farc;
        MPI_Type_vector(inner_cnt, blklen, stride, MPI_DOUBLE, &newtype_farc);
        MPI_Type_commit(&newtype_farc);
        HRT_GET_TIMESTAMP(stop);
        HRT_GET_ELAPSED_TICKS(start, stop, &interposer_type_create);

        for (int i=0; i<inner_runs; i++) {
            HRT_GET_TIMESTAMP(start);
            if (rank == 0) {
                farc::DDT_Pack(farc_inbuf, farc_outbuf, t2, outer_cnt);
                PMPI_Send(farc_outbuf, t2_size, MPI_BYTE, 1, 0, MPI_COMM_WORLD);
            }
            else {
                PMPI_Recv(farc_outbuf, t2_size, MPI_BYTE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            HRT_GET_TIMESTAMP(stop);
            HRT_GET_ELAPSED_TICKS(start, stop, &farc_send);
            MPI_Barrier(MPI_COMM_WORLD);
        }

        for (int i=0; i<inner_runs; i++) {
            HRT_GET_TIMESTAMP(start);
            if (rank == 0) {
                MPI_Send(interposer_inbuf, outer_cnt, newtype_farc, 1, 0, MPI_COMM_WORLD);
            }
            else {
                MPI_Recv(interposer_outbuf, t2_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            HRT_GET_TIMESTAMP(stop);
            HRT_GET_ELAPSED_TICKS(start, stop, &interposer_send);
            MPI_Barrier(MPI_COMM_WORLD);
        }

        for (int i=0; i<inner_runs; i++) {
            HRT_GET_TIMESTAMP(start);
            if (rank == 0) {
                PMPI_Send(mpi_inbuf, outer_cnt, newtype, 1, 0, MPI_COMM_WORLD);
            }
            else {
                PMPI_Recv(mpi_outbuf, t2_size, MPI_BYTE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            HRT_GET_TIMESTAMP(stop);
            HRT_GET_ELAPSED_TICKS(start, stop, &mpi_send);
            MPI_Barrier(MPI_COMM_WORLD);
        }
	      
        for (int i=0; i<inner_runs; i++) {
            HRT_GET_TIMESTAMP(start);
            if (rank == 0) {
                for(j=0; j < inner_cnt; ++j) {
                    cpp_outbuf[j*5    ] = cpp_inbuf[j*stride*5];
                    cpp_outbuf[j*5 + 1] = cpp_inbuf[j*stride*5 + 1];
                    cpp_outbuf[j*5 + 2] = cpp_inbuf[j*stride*5 + 2];
                    cpp_outbuf[j*5 + 3] = cpp_inbuf[j*stride*5 + 3];
                    cpp_outbuf[j*5 + 4] = cpp_inbuf[j*stride*5 + 4];
                }
                PMPI_Send(cpp_outbuf, t2_size, MPI_BYTE, 1, 0, MPI_COMM_WORLD);
            }
            else {
                PMPI_Recv(cpp_outbuf, t2_size, MPI_BYTE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            HRT_GET_TIMESTAMP(stop);
            HRT_GET_ELAPSED_TICKS(start, stop, &cpp_send)
            MPI_Barrier(MPI_COMM_WORLD);
        }

        if (rank == 0) {
            static int firstline=1;
            if (firstline) printf("%10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n", "size", "mpi_create", "psr_create", "farc_create", "cpp_send", "mpi_send", "psr_send", "farc_send", "blklen", "stride", "count", "send_count");
            firstline=0;
            printf("%10i %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10i %10i %10i %10i\n", data_size, HRT_GET_USEC(mpi_type_create), HRT_GET_USEC(interposer_type_create), HRT_GET_USEC(farc_type_create), HRT_GET_USEC(cpp_send), HRT_GET_USEC(mpi_send), HRT_GET_USEC(interposer_send), HRT_GET_USEC(farc_send), blklen, stride, inner_cnt, outer_cnt);
        }

        PMPI_Type_free(&newtype);
        MPI_Type_free(&newtype_farc);
        farc::DDT_Free(t1);
        farc::DDT_Free(t2);

	  }
	
    if (rank == 1) {
        int res = compare_buffers(buffer_size, &mpi_inbuf, &farc_inbuf, &mpi_outbuf, &farc_outbuf);
        test_result(res);
    }

//    free_buffers(&mpi_inbuf, &farc_inbuf, &mpi_outbuf, &farc_outbuf);

}

int main(int argc, char** argv) {

   if (argc < 14) {
        fprintf(stderr, "%s [num-runs] [blklen_start] [blklen_end] [blklen_inc] [stride_start] [stride_end] [stride_inc] [inner_cnt_start] [inner_cnt_end] [inner_cnt_inc] [outer_cnt_start] [outer_cnt_end] [outer_cnt_inc]\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    MPI_Init(&argc, &argv);
    HRT_INIT(1, g_timerfreq);
 
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

