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

    int buffer_size = 8 * sizeof(double);
    int sendnum = 4;

    char* mpi_inbuf;
    char* mpi_outbuf;
    char* pmpi_inbuf;
    char* pmpi_outbuf;

    test_start("send/irecv bug test");
    init_buffers(buffer_size, &mpi_inbuf, &pmpi_inbuf, &mpi_outbuf, &pmpi_outbuf);

    MPI_Datatype hvec1_t;
    MPI_Datatype hvec_t;
    MPI_Type_hvector(2, 1, 16, MPI_DOUBLE, &hvec1_t);
    MPI_Type_hvector(2, 1, 32, hvec1_t, &hvec_t);
    MPI_Type_commit(&hvec_t);

    MPI_Datatype pmpi_hvec1_t;
    MPI_Datatype pmpi_hvec_t;
    PMPI_Type_hvector(2, 1, 16, MPI_DOUBLE, &pmpi_hvec1_t);
    PMPI_Type_hvector(2, 1, 32, pmpi_hvec1_t, &pmpi_hvec_t);
    PMPI_Type_commit(&pmpi_hvec_t);



    if (rank % 2 == 0) {
        MPI_Send(mpi_inbuf, sendnum, MPI_DOUBLE, peer, 0, MPI_COMM_WORLD);
        MPI_Recv(mpi_outbuf, 1, hvec_t, peer, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        PMPI_Send(pmpi_inbuf, sendnum, MPI_DOUBLE, peer, 0, MPI_COMM_WORLD);
        PMPI_Recv(pmpi_outbuf, 1, pmpi_hvec_t, peer, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else {
        MPI_Recv(mpi_outbuf, 1, hvec_t, peer, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(mpi_inbuf, sendnum, MPI_DOUBLE, peer, 0, MPI_COMM_WORLD);

        PMPI_Recv(pmpi_outbuf, 1, pmpi_hvec_t, peer, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        PMPI_Send(pmpi_inbuf, sendnum, MPI_DOUBLE, peer, 0, MPI_COMM_WORLD);
    }

    int res = compare_buffers(buffer_size, &mpi_inbuf, &pmpi_inbuf, &mpi_outbuf, &pmpi_outbuf);
    test_result(res);

    free_buffers(&mpi_inbuf, &pmpi_inbuf, &mpi_outbuf, &pmpi_outbuf);
    MPI_Finalize();
}

