// Copyright 2013 Timo Schneider and Fredrik Berg Kjolstad
//
// This file is part of the libpack packing library.
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT in the top level directory for details.

#include <string>
#include <mpi.h>

#include "test.hpp"
#include "../ddt_jit.hpp"

int main(int argc, char** argv) {

    char* mpi_inbuf;
    char* mpi_outbuf;
    char* farc_inbuf;
    char* farc_outbuf;

    MPI_Init(&argc, &argv);

    test_start("pack(2, resized(lb=0, extent=2, [ctg(2)[MPI_INT]]))");
    init_buffers(20*sizeof(int), &mpi_inbuf, &farc_inbuf, &mpi_outbuf, &farc_outbuf);

    farc::DDT_Init();
    farc::Datatype* t1 = new farc::PrimitiveDatatype(farc::PrimitiveDatatype::INT);
    farc::Datatype* t2 = new farc::ContiguousDatatype(2, t1);
    farc::Datatype* t3 = new farc::ResizedDatatype(t2, 0, 2);


    farc::DDT_Commit(t3);
    farc::DDT_Pack(farc_inbuf, farc_outbuf, t3, 2);

    int position = 0;
    MPI_Datatype mpiresized, mpicont;
    MPI_Type_contiguous(2, MPI_INT, &mpicont);
    MPI_Type_create_resized(mpicont, 0, 2, &mpiresized);
    MPI_Type_commit(&mpiresized);
    
    int res = compare_ddt_info(mpiresized, t3);
    
    MPI_Pack(mpi_inbuf, 2, mpiresized, mpi_outbuf, 20*sizeof(int), &position, MPI_COMM_WORLD);

    res += compare_buffers(20*sizeof(int), &mpi_inbuf, &farc_inbuf, &mpi_outbuf, &farc_outbuf);
    // inspect_buffers(20*sizeof(int), &mpi_inbuf, &farc_inbuf, &mpi_outbuf, &farc_outbuf);
    free_buffers(&mpi_inbuf, &farc_inbuf, &mpi_outbuf, &farc_outbuf);
    test_result(res);

    MPI_Finalize();

    return 0;

}

