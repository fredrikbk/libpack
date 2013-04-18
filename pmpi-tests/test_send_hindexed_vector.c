// Copyright 2013 Timo Schneider and Fredrik Berg Kjolstad
//
// This file is part of the libpack packing library.
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT in the top level directory for details.

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

    test_start("send/recv (2, hindexed[{(1*vector[[int], count=2, blklen=1, stride=2], offset=0), (1*vector[[int], count=2, blklen=1, stride=2], offset=40)}])");
    init_buffers(20*sizeof(int), &mpi_inbuf, &pmpi_inbuf, &mpi_outbuf, &pmpi_outbuf);

    int blocklen[2] = {1, 1};
    MPI_Aint disp[2] = {0, 40};

    MPI_Datatype subtype;
    MPI_Type_vector(2, 1, 2, MPI_INT, &subtype);
    MPI_Datatype mpitype; 
    MPI_Type_create_hindexed(2, blocklen, disp, subtype, &mpitype);
    MPI_Type_commit(&mpitype);

    MPI_Datatype pmpi_subtype;
    PMPI_Type_vector(2, 1, 2, MPI_INT, &pmpi_subtype);
    MPI_Datatype pmpitype; 
    PMPI_Type_create_hindexed(2, blocklen, disp, pmpi_subtype, &pmpitype);
    PMPI_Type_commit(&pmpitype);

    if (rank % 2 == 0) {
        MPI_Send(mpi_inbuf, 1, mpitype, peer, 0, MPI_COMM_WORLD);
        MPI_Recv(mpi_outbuf, 1, mpitype, peer, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);       

        PMPI_Send(pmpi_inbuf, 1, pmpitype, peer, 0, MPI_COMM_WORLD);
        PMPI_Recv(pmpi_outbuf, 1, pmpitype, peer, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);       
    }
    else {
        MPI_Recv(mpi_outbuf, 1, mpitype, peer, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);       
        MPI_Send(mpi_inbuf, 1, mpitype, peer, 0, MPI_COMM_WORLD);

        PMPI_Recv(pmpi_outbuf, 1, pmpitype, peer, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);       
        PMPI_Send(pmpi_inbuf, 1, pmpitype, peer, 0, MPI_COMM_WORLD);
    }

    int res = compare_buffers(20*sizeof(int), &mpi_inbuf, &pmpi_inbuf, &mpi_outbuf, &pmpi_outbuf);
    free_buffers(&mpi_inbuf, &pmpi_inbuf, &mpi_outbuf, &pmpi_outbuf);
    test_result(res);

    MPI_Finalize();

}

