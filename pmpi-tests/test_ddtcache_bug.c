#include "test.hpp"
#include <mpi.h>

// this number should be bigger than the interposers cache size
#define NUM_DDTS_TO_TEST 1000

int main(int argc, char** argv) {

    MPI_Init(&argc, &argv);
    int rank, peer, commsize;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    if (rank % 2) peer = rank - 1 % commsize;
    else peer = rank + 1 % commsize;
    
    if (commsize % 2 != 0) {
        fprintf(stderr, "Use even number of processes.\n");
        exit(EXIT_FAILURE);
    }

    int buffer_size = 8 * sizeof(double);
    int sendnum = 4;

    char* mpi_inbuf;
    char* mpi_outbuf;
    char* pmpi_inbuf;
    char* pmpi_outbuf;

    test_start("ddt cache bug");
    init_buffers(buffer_size, &mpi_inbuf, &pmpi_inbuf, &mpi_outbuf, &pmpi_outbuf);

    MPI_Datatype* mpiddts;
    MPI_Datatype* pmpiddts;
    mpiddts = (MPI_Datatype*) malloc(NUM_DDTS_TO_TEST * sizeof(MPI_Datatype));
    pmpiddts = (MPI_Datatype*) malloc(NUM_DDTS_TO_TEST * sizeof(MPI_Datatype));

    int i=0;
    int res = 0;
    for (i=0; i<NUM_DDTS_TO_TEST; i++) {
        MPI_Type_contiguous(sendnum, MPI_DOUBLE, &(mpiddts[i]));
        MPI_Type_commit(&(mpiddts[i]));
        PMPI_Type_contiguous(sendnum, MPI_DOUBLE, &(pmpiddts[i]));
        PMPI_Type_commit(&(pmpiddts[i]));
 

        if (rank % 2 == 0) {
            MPI_Send(mpi_inbuf, sendnum, MPI_DOUBLE, peer, 0, MPI_COMM_WORLD);
            MPI_Recv(mpi_outbuf, 1, mpiddts[i], peer, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            PMPI_Send(pmpi_inbuf, sendnum, MPI_DOUBLE, peer, 0, MPI_COMM_WORLD);
            PMPI_Recv(pmpi_outbuf, 1, pmpiddts[i], peer, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        else {
            MPI_Recv(mpi_outbuf, 1, mpiddts[i], peer, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(mpi_inbuf, sendnum, MPI_DOUBLE, peer, 0, MPI_COMM_WORLD);

            PMPI_Recv(pmpi_outbuf, 1, pmpiddts[i], peer, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            PMPI_Send(pmpi_inbuf, sendnum, MPI_DOUBLE, peer, 0, MPI_COMM_WORLD);
        }

        res += compare_buffers(buffer_size, &mpi_inbuf, &pmpi_inbuf, &mpi_outbuf, &pmpi_outbuf);
    }
    test_result(res);

    free_buffers(&mpi_inbuf, &pmpi_inbuf, &mpi_outbuf, &pmpi_outbuf);

    MPI_Finalize();
}

