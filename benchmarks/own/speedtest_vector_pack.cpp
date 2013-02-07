#include <string>

#include "../../ddt_jit.hpp"
#include "../../tests/test.hpp"
#include "../../copy_benchmark/hrtimer/hrtimer.h"

unsigned long long g_timerfreq;

void benchmark_vector(int blklen, int stride, int inner_cnt, int outer_cnt, int inner_runs, int outer_runs) {

    char* mpi_inbuf;
    char* mpi_outbuf;
    char* farc_inbuf;
    char* farc_outbuf;
    HRT_TIMESTAMP_T start, stop;
    uint64_t mpi_type_create, farc_type_create, mpi_pack, farc_pack;

    int data_size = sizeof(int)*inner_cnt * outer_cnt * blklen;
    int buffer_size = sizeof(int)*((inner_cnt-1)*stride+blklen) * outer_cnt;

    init_buffers(buffer_size, &mpi_inbuf, &farc_inbuf, &mpi_outbuf, &farc_outbuf);

    FARC_DDT_Init();

    for (int o=0; o<outer_runs; o++) {
	    HRT_GET_TIMESTAMP(start);
	    FARC_Datatype* t1 = new FARC_PrimitiveDatatype(MPI_INT);
	    FARC_Datatype* t2 = new FARC_VectorDatatype(t1, inner_cnt, blklen, stride);
	    int ddt_handle = FARC_DDT_Commit(t2);
	    HRT_GET_TIMESTAMP(stop);
	    HRT_GET_ELAPSED_TICKS(start, stop, &farc_type_create);
	
	    HRT_GET_TIMESTAMP(start);
	    MPI_Datatype newtype;
	    MPI_Type_vector(inner_cnt, blklen, stride, MPI_INT, &newtype);
	    MPI_Type_commit(&newtype);
	    HRT_GET_TIMESTAMP(stop);
	    HRT_GET_ELAPSED_TICKS(start, stop, &mpi_type_create);
	
	    for (int i=0; i<inner_runs; i++) {
	
	        HRT_GET_TIMESTAMP(start);
	        FARC_DDT_Pack(farc_inbuf, farc_outbuf, ddt_handle, outer_cnt);
	        HRT_GET_TIMESTAMP(stop);
	        HRT_GET_ELAPSED_TICKS(start, stop, &farc_pack);

	        HRT_GET_TIMESTAMP(start);
	        int position = 0;
	        MPI_Pack(mpi_inbuf, outer_cnt, newtype, mpi_outbuf, buffer_size*sizeof(int), &position, MPI_COMM_WORLD);
	        HRT_GET_TIMESTAMP(stop);
	        HRT_GET_ELAPSED_TICKS(start, stop, &mpi_pack);

    	    static int firstline=1;
    	    if (firstline) printf("%10s %10s %10s %10s %10s %10s %10s %10s %10s\n", "size", "mpi_create", "farc_create", "mpi_pack", "farc_pack", "blklen", "stride", "count", "pack_count");
            firstline=0;
    	    printf("%10i %10.3lf %10.3lf %10.3lf %10.3lf %10i %10i %10i %10i\n", data_size, HRT_GET_USEC(mpi_type_create), HRT_GET_USEC(farc_type_create), HRT_GET_USEC(mpi_pack), HRT_GET_USEC(farc_pack), blklen, stride, inner_cnt, outer_cnt);
	    } 
	
    }	

	free_buffers(&mpi_inbuf, &farc_inbuf, &mpi_outbuf, &farc_outbuf);

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

