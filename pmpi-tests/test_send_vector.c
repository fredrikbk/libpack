#include "test.hpp"

#include <mpi.h>

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

    char* mpi_inbuf;
    char* mpi_outbuf;
    char* pmpi_inbuf;
    char* pmpi_outbuf;

    test_start("send/recv (2, vector[[int], count=2, blklen=3, stride=5])");
    init_buffers(20*sizeof(int), &mpi_inbuf, &pmpi_inbuf, &mpi_outbuf, &pmpi_outbuf);

    MPI_Datatype vector_ddt;
    MPI_Type_vector(2, 3, 5, MPI_INT, &vector_ddt);
    MPI_Type_commit(&vector_ddt);

    if (rank % 2 == 0) {
        MPI_Send(mpi_inbuf, 2, vector_ddt, peer, 0, MPI_COMM_WORLD);
        MPI_Recv(mpi_outbuf, 2, vector_ddt, peer, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);       

        PMPI_Send(pmpi_inbuf, 2, vector_ddt, peer, 0, MPI_COMM_WORLD);
        PMPI_Recv(pmpi_outbuf, 2, vector_ddt, peer, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);       
    }
    else {
        MPI_Recv(mpi_outbuf, 2, vector_ddt, peer, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);       
        MPI_Send(mpi_inbuf, 2, vector_ddt, peer, 0, MPI_COMM_WORLD);

        PMPI_Recv(pmpi_outbuf, 2, vector_ddt, peer, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);       
        PMPI_Send(pmpi_inbuf, 2, vector_ddt, peer, 0, MPI_COMM_WORLD);
    }

    int res = compare_buffers(20*sizeof(int), &mpi_inbuf, &pmpi_inbuf, &mpi_outbuf, &pmpi_outbuf);
    free_buffers(&mpi_inbuf, &pmpi_inbuf, &mpi_outbuf, &pmpi_outbuf);
    test_result(res);

    MPI_Finalize();

}

