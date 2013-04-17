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

    test_start("send/recv (2, hindexed[{(1*MPI_INT, offset=4), (3*MPI_INT, offset=16), (2*MPI_INT, offset=32)}])");
    init_buffers(20*sizeof(int), &mpi_inbuf, &pmpi_inbuf, &mpi_outbuf, &pmpi_outbuf);

    int blocklen[3] = {1, 3, 2};
    MPI_Aint disp[3] = {4, 16, 32};

    MPI_Datatype mpitype; 
    MPI_Type_create_hindexed(3, blocklen, disp, MPI_INT, &mpitype);
    MPI_Type_commit(&mpitype);

    MPI_Datatype pmpitype; 
    PMPI_Type_create_hindexed(3, blocklen, disp, MPI_INT, &pmpitype);
    PMPI_Type_commit(&pmpitype);

    if (rank % 2 == 0) {
        MPI_Send(mpi_inbuf, 2, mpitype, peer, 0, MPI_COMM_WORLD);
        MPI_Recv(mpi_outbuf, 2, mpitype, peer, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);       

        PMPI_Send(pmpi_inbuf, 2, pmpitype, peer, 0, MPI_COMM_WORLD);
        PMPI_Recv(pmpi_outbuf, 2, pmpitype, peer, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);       
    }
    else {
        MPI_Recv(mpi_outbuf, 2, mpitype, peer, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);       
        MPI_Send(mpi_inbuf, 2, mpitype, peer, 0, MPI_COMM_WORLD);

        PMPI_Recv(pmpi_outbuf, 2, pmpitype, peer, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);       
        PMPI_Send(pmpi_inbuf, 2, pmpitype, peer, 0, MPI_COMM_WORLD);
    }

    int res = compare_buffers(20*sizeof(int), &mpi_inbuf, &pmpi_inbuf, &mpi_outbuf, &pmpi_outbuf);
    free_buffers(&mpi_inbuf, &pmpi_inbuf, &mpi_outbuf, &pmpi_outbuf);
    test_result(res);

    MPI_Finalize();

}

