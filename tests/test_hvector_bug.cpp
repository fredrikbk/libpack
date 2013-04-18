// Copyright 2013 Timo Schneider and Fredrik Berg Kjolstad
//
// This file is part of the libpack packing library.
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT in the top level directory for details.

#include <string>
#include <mpi.h>

#include "../ddt_jit.hpp"
#include "test.hpp"

int main(int argc, char** argv) {

    double* mpi_inbuf;
    double* mpi_outbuf;
    double* farc_inbuf;
    double* farc_outbuf;

    MPI_Init(&argc, &argv);

    int buffer_size = 414720*sizeof(double);

    int count    = 64;
    int blocklen = 1;
    int stride   = 51840;

    test_start("hvector bug");

    mpi_inbuf = (double*) malloc(buffer_size);
    mpi_outbuf = (double*) malloc(buffer_size);
    farc_inbuf = (double*) malloc(buffer_size);
    farc_outbuf = (double*) malloc(buffer_size);

    for (size_t i=0; i<buffer_size/sizeof(double); i++) {
        mpi_inbuf[i] = i+1;
        farc_inbuf[i] = i+1;
        mpi_outbuf[i] = 0;
        farc_outbuf[i] = 0;
    }

    farc::DDT_Init();
    farc::Datatype* t1 = new farc::PrimitiveDatatype(farc::PrimitiveDatatype::DOUBLE);
    farc::Datatype* t2 = new farc::HVectorDatatype(count, blocklen, stride, t1);
    farc::DDT_Commit(t2);
    farc::DDT_Pack((char*)farc_inbuf, (char*)farc_outbuf, t2, 1);

    MPI_Datatype newtype;
    MPI_Type_hvector(count, blocklen, stride, MPI_DOUBLE, &newtype);
    MPI_Type_commit(&newtype);
    int position = 0;
    MPI_Pack(mpi_inbuf, 1, newtype, mpi_outbuf, buffer_size, &position, MPI_COMM_WORLD);

    int res = compare_buffers(buffer_size, (char**)&mpi_inbuf, (char**)&farc_inbuf, (char**)&mpi_outbuf, (char**)&farc_outbuf);
    test_result(res);

    /*
    printf("\n");
    for (int i=0; i<buffer_size/sizeof(double); i++) {
        printf("mpi_inbuf[%i] = %lf farc_inbuf[%i] = %lf mpi_outbuf[%i] = %lf farc_outbuf[%i] = %lf\n", i, mpi_inbuf[i], i, farc_inbuf[i], i, mpi_outbuf[i], i, farc_outbuf[i]);
    }
    */
    free(mpi_inbuf);
    free(mpi_outbuf);
    free(farc_inbuf);
    free(farc_outbuf);


    MPI_Finalize();

    return 0;

}

