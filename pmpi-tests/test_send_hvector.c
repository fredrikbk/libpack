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

    test_start("send/recv (2, hvector[[int], count=2, blklen=3, stride=5*sizeof(int)])");
    init_buffers(20*sizeof(int), &mpi_inbuf, &pmpi_inbuf, &mpi_outbuf, &pmpi_outbuf);

    MPI_Datatype hvector_ddt;
    MPI_Type_hvector(2, 3, 5*sizeof(int), MPI_INT, &hvector_ddt);
    MPI_Type_commit(&hvector_ddt);

    MPI_Datatype pmpi_hvector_ddt;
    PMPI_Type_hvector(2, 3, 5*sizeof(int), MPI_INT, &pmpi_hvector_ddt);
    PMPI_Type_commit(&pmpi_hvector_ddt);

    if (rank % 2 == 0) {
        MPI_Send(mpi_inbuf, 2, hvector_ddt, peer, 0, MPI_COMM_WORLD);
        MPI_Recv(mpi_outbuf, 2, hvector_ddt, peer, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);       

        PMPI_Send(pmpi_inbuf, 2, pmpi_hvector_ddt, peer, 0, MPI_COMM_WORLD);
        PMPI_Recv(pmpi_outbuf, 2, pmpi_hvector_ddt, peer, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);       
    }
    else {
        MPI_Recv(mpi_outbuf, 2, hvector_ddt, peer, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);       
        MPI_Send(mpi_inbuf, 2, hvector_ddt, peer, 0, MPI_COMM_WORLD);

        PMPI_Recv(pmpi_outbuf, 2, pmpi_hvector_ddt, peer, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);       
        PMPI_Send(pmpi_inbuf, 2, pmpi_hvector_ddt, peer, 0, MPI_COMM_WORLD);
    }

    int res = compare_buffers(20*sizeof(int), &mpi_inbuf, &pmpi_inbuf, &mpi_outbuf, &pmpi_outbuf);
    free_buffers(&mpi_inbuf, &pmpi_inbuf, &mpi_outbuf, &pmpi_outbuf);
    test_result(res);

    MPI_Finalize();

}

