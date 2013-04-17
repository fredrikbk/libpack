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

    test_start("isend/irecv + wait (2, vector[[int], count=2, blklen=3, stride=5])");
    init_buffers(20*sizeof(int), &mpi_inbuf, &pmpi_inbuf, &mpi_outbuf, &pmpi_outbuf);

    MPI_Datatype vector_ddt;
    MPI_Type_vector(2, 3, 5, MPI_INT, &vector_ddt);
    MPI_Type_commit(&vector_ddt);

    MPI_Datatype pmpi_vector_ddt;
    PMPI_Type_vector(2, 3, 5, MPI_INT, &pmpi_vector_ddt);
    PMPI_Type_commit(&pmpi_vector_ddt);

    MPI_Request requests_mpi[2];
    MPI_Request requests_pmpi[2];
    MPI_Status statuses_mpi[2]; 
    MPI_Status statuses_pmpi[2];

    if (rank % 2 == 0) {
        MPI_Isend(mpi_inbuf, 2, vector_ddt, peer, 0, MPI_COMM_WORLD, &(requests_mpi[0]));
        MPI_Irecv(mpi_outbuf, 2, vector_ddt, peer, 0, MPI_COMM_WORLD, &(requests_mpi[1]));

        PMPI_Isend(pmpi_inbuf, 2, pmpi_vector_ddt, peer, 0, MPI_COMM_WORLD, &(requests_pmpi[0]));
        PMPI_Irecv(pmpi_outbuf, 2, pmpi_vector_ddt, peer, 0, MPI_COMM_WORLD, &(requests_pmpi[1]));       
    }
    else {
        MPI_Irecv(mpi_outbuf, 2, vector_ddt, peer, 0, MPI_COMM_WORLD, &(requests_mpi[0]));       
        MPI_Isend(mpi_inbuf, 2, vector_ddt, peer, 0, MPI_COMM_WORLD, &(requests_mpi[1]));

        PMPI_Irecv(pmpi_outbuf, 2, pmpi_vector_ddt, peer, 0, MPI_COMM_WORLD, &(requests_pmpi[0]));       
        PMPI_Isend(pmpi_inbuf, 2, pmpi_vector_ddt, peer, 0, MPI_COMM_WORLD, &(requests_pmpi[1]));
    }

    MPI_Wait(&(requests_mpi[0]),  &(statuses_mpi[0]));
    MPI_Wait(&(requests_mpi[1]),  &(statuses_mpi[1]));
    MPI_Wait(&(requests_pmpi[0]), &(statuses_pmpi[0]));
    MPI_Wait(&(requests_pmpi[1]), &(statuses_pmpi[1]));

    int res = compare_buffers(20*sizeof(int), &mpi_inbuf, &pmpi_inbuf, &mpi_outbuf, &pmpi_outbuf);
    free_buffers(&mpi_inbuf, &pmpi_inbuf, &mpi_outbuf, &pmpi_outbuf);
    test_result(res);

    MPI_Finalize();

}

