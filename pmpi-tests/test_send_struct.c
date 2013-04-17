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

    struct Partstruct { 
        int foo;
        int moo[6];
        char bar[8];
    };

    struct Partstruct* particle;
    particle = (struct Partstruct*) mpi_inbuf;

    test_start("send/recv (2, struct[{1*MPI_INT, 6*MPI_INT, 8*MPI_CHAR}])");
    init_buffers(20*sizeof(int), &mpi_inbuf, &pmpi_inbuf, &mpi_outbuf, &pmpi_outbuf);

    MPI_Datatype type[3] = {MPI_INT, MPI_INT, MPI_CHAR};
    int blocklen[3] = {1, 6, 8};
    MPI_Aint disp[3];
    int base;

    // calc displacements
    MPI_Address(&(particle[0].foo), disp+0); 
    MPI_Address(&(particle[0].moo), disp+1); 
    MPI_Address(&(particle[0].bar), disp+2); 
    base = disp[0]; 
    int i;
    for (i=0; i <3; i++) disp[i] -= base;

    MPI_Datatype mpitype; 
    MPI_Type_struct(3, blocklen, disp, type, &mpitype);
    MPI_Type_commit(&mpitype);

    MPI_Datatype pmpitype; 
    PMPI_Type_struct(3, blocklen, disp, type, &pmpitype);
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

