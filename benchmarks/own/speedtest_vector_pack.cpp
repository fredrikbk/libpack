#include <string>

#include <ddt_jit.hpp>
#include "../../tests/test.hpp"
#include "../../copy_benchmark/hrtimer/hrtimer.h"


unsigned long long g_timerfreq;

void init_in_and_out_buffer(size_t buffer_size, char** inbuf, char** outbuf) {

	*inbuf = (char*) malloc(buffer_size);
    *outbuf = (char*) malloc(buffer_size);

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
    char* pmpi_inbuf;
    char* pmpi_outbuf;
    char* farc_inbuf;
    char* farc_outbuf;
    double* cpp_inbuf;
    double* cpp_outbuf;
    int jend, j;

    HRT_TIMESTAMP_T start, stop;
    uint64_t mpi_type_create, pmpi_type_create, farc_type_create, mpi_pack, pmpi_pack, farc_pack, cpp_pack;

    int data_size = sizeof(double)*inner_cnt * outer_cnt * blklen;
    int buffer_size = sizeof(double)*((inner_cnt-1)*stride+blklen) * outer_cnt;

	init_in_and_out_buffer(buffer_size, &mpi_inbuf, &mpi_outbuf);
	init_in_and_out_buffer(buffer_size, &pmpi_inbuf, &pmpi_outbuf);
	init_in_and_out_buffer(buffer_size, &farc_inbuf, &farc_outbuf);
	init_in_and_out_buffer(buffer_size, reinterpret_cast<char**>(&cpp_inbuf),  reinterpret_cast<char**>(&cpp_outbuf));

    for (int o=0; o<outer_runs; o++) {
        HRT_GET_TIMESTAMP(start);
        FARC_Datatype* t1 = new FARC_PrimitiveDatatype(MPI_DOUBLE);
        FARC_Datatype* t2 = new FARC_VectorDatatype(t1, inner_cnt, blklen, stride);
        FARC_DDT_Commit(t2);
        HRT_GET_TIMESTAMP(stop);
        HRT_GET_ELAPSED_TICKS(start, stop, &farc_type_create);

        HRT_GET_TIMESTAMP(start);
        MPI_Datatype newtype_mpi;
        MPI_Type_vector(inner_cnt, blklen, stride, MPI_DOUBLE, &newtype_mpi);
        MPI_Type_commit(&newtype_mpi);
        HRT_GET_TIMESTAMP(stop);
        HRT_GET_ELAPSED_TICKS(start, stop, &mpi_type_create);

        HRT_GET_TIMESTAMP(start);
        MPI_Datatype newtype_pmpi;
        PMPI_Type_vector(inner_cnt, blklen, stride, MPI_DOUBLE, &newtype_pmpi);
        PMPI_Type_commit(&newtype_pmpi);
        HRT_GET_TIMESTAMP(stop);
        HRT_GET_ELAPSED_TICKS(start, stop, &pmpi_type_create);

        for (int i=0; i<inner_runs; i++) {
            int position = 0;
            HRT_GET_TIMESTAMP(start);
            MPI_Pack(mpi_inbuf, outer_cnt, newtype_mpi, mpi_outbuf, buffer_size*sizeof(int), &position, MPI_COMM_WORLD);
            HRT_GET_TIMESTAMP(stop);
            HRT_GET_ELAPSED_TICKS(start, stop, &mpi_pack);

            position = 0;
            HRT_GET_TIMESTAMP(start);
            PMPI_Pack(pmpi_inbuf, outer_cnt, newtype_pmpi, mpi_outbuf, buffer_size*sizeof(int), &position, MPI_COMM_WORLD);
            HRT_GET_TIMESTAMP(stop);
            HRT_GET_ELAPSED_TICKS(start, stop, &pmpi_pack);

            HRT_GET_TIMESTAMP(start);
            for (int i=0; i < outer_cnt; i++) {
                for(int j=0; j < inner_cnt; j++) {
                    cpp_outbuf[i*inner_cnt*5 + j*5    ] = cpp_inbuf[i*((inner_cnt-1)*stride + 5) + j*stride*5];
                    cpp_outbuf[i*inner_cnt*5 + j*5 + 1] = cpp_inbuf[i*((inner_cnt-1)*stride + 5) + j*stride*5 + 1];
                    cpp_outbuf[i*inner_cnt*5 + j*5 + 2] = cpp_inbuf[i*((inner_cnt-1)*stride + 5) + j*stride*5 + 2];
                    cpp_outbuf[i*inner_cnt*5 + j*5 + 3] = cpp_inbuf[i*((inner_cnt-1)*stride + 5) + j*stride*5 + 3];
                    cpp_outbuf[i*inner_cnt*5 + j*5 + 4] = cpp_inbuf[i*((inner_cnt-1)*stride + 5) + j*stride*5 + 4];
                }
            }
            HRT_GET_TIMESTAMP(stop);
            HRT_GET_ELAPSED_TICKS(start, stop, &cpp_pack);

            HRT_GET_TIMESTAMP(start);
            FARC_DDT_Pack(farc_inbuf, farc_outbuf, t2, outer_cnt);
            HRT_GET_TIMESTAMP(stop);
            HRT_GET_ELAPSED_TICKS(start, stop, &farc_pack);

            static int firstline=1;
            if (firstline) printf("%10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n", "size", "mpi_create", "pmpi_create", "farc_create", "cpp_pack", "mpi_pack", "pmpi_pack", "farc_pack", "blklen", "stride", "count", "pack_count");
            firstline=0;
            printf("%10i %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10i %10i %10i %10i\n", data_size, HRT_GET_USEC(mpi_type_create), HRT_GET_USEC(pmpi_type_create), HRT_GET_USEC(farc_type_create), HRT_GET_USEC(cpp_pack), HRT_GET_USEC(mpi_pack), HRT_GET_USEC(pmpi_pack), HRT_GET_USEC(farc_pack), blklen, stride, inner_cnt, outer_cnt);

        } 
        MPI_Type_free(&newtype_mpi);
        MPI_Type_free(&newtype_pmpi);
        FARC_DDT_Free(t1);
        FARC_DDT_Free(t2);

	  }
	
//    int res = compare_buffers(buffer_size, &mpi_inbuf, &farc_inbuf, &mpi_outbuf, &farc_outbuf);
//    test_result(res);

	free_in_and_out_buffer(buffer_size, &mpi_inbuf, &mpi_outbuf);
	free_in_and_out_buffer(buffer_size, &pmpi_inbuf, &pmpi_outbuf);
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
    FARC_DDT_Init();
 
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

