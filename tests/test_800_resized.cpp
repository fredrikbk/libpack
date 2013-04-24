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

    char* mpi_inbuf;
    char* mpi_outbuf;
    char* farc_inbuf;
    char* farc_outbuf;

    MPI_Init(&argc, &argv);

    test_start("pack(2, resize(lb=0, ext=2, [vector([int], count=3, blklen=1, stride=3)]))");
    init_buffers(20*sizeof(int), &mpi_inbuf, &farc_inbuf, &mpi_outbuf, &farc_outbuf);

    farc::DDT_Init();
    farc::Datatype* t1 = new farc::PrimitiveDatatype(farc::PrimitiveDatatype::INT);
    farc::Datatype* t2 = new farc::VectorDatatype(3, 1, 3, t1);
    farc::Datatype* t3 = new farc::ResizedDatatype(t2, 0, 2);
    farc::DDT_Commit(t3);
    farc::DDT_Pack(farc_inbuf, farc_outbuf, t3, 2);

    MPI_Datatype vectype, newtype;
    MPI_Type_vector(3, 1, 3, MPI_INT, &vectype);
    MPI_Type_create_resized(vectype, 0, 2, &newtype);
    MPI_Type_commit(&newtype);
    int position = 0;
    MPI_Pack(mpi_inbuf, 2, newtype, mpi_outbuf, 20*sizeof(int), &position, MPI_COMM_WORLD);

    int res = compare_ddt_info(newtype, t3);
    res += compare_buffers(20*sizeof(int), &mpi_inbuf, &farc_inbuf, &mpi_outbuf, &farc_outbuf);
    // inspect_buffers(20*sizeof(int), &mpi_inbuf, &farc_inbuf, &mpi_outbuf, &farc_outbuf);
    
    free_buffers(&mpi_inbuf, &farc_inbuf, &mpi_outbuf, &farc_outbuf);
    test_result(res);

    MPI_Finalize();

    return 0;

}

