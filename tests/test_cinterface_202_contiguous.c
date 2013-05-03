// Copyright 2013 Timo Schneider and Fredrik Berg Kjolstad
//
// This file is part of the libpack packing library.
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT in the top level directory for details.

#include <mpi.h>

#include "test.h"
#include "../pack.h"

int main(int argc, char** argv) {

    char* mpi_inbuf;
    char* mpi_outbuf;
    char* farc_inbuf;
    char* farc_outbuf;

    MPI_Init(&argc, &argv);
    LPK_Init();

    test_start("pack(2, contiguous[[int], count=2]) [using c interface]");
    init_buffers(20*sizeof(int), &mpi_inbuf, &farc_inbuf, &mpi_outbuf, &farc_outbuf);

//  farc::Datatype* t1 = new farc::PrimitiveDatatype(farc::PrimitiveDatatype::INT);
    LPK_Datatype t1;
    LPK_Primitive(LPK_INT, &t1);        

//  farc::Datatype* t2 = new farc::ContiguousDatatype(0, t1);
    LPK_Datatype t2;
    LPK_Contiguous(2, t1, &t2);

//  farc::DDT_Commit(t2);
    LPK_Compile(t2);

//  farc::DDT_Pack(farc_inbuf, farc_outbuf, t2, 2);
    LPK_Pack(farc_inbuf, 2, t2, farc_outbuf);

    int position = 0;
    MPI_Datatype new_contig;
    MPI_Type_contiguous(2, MPI_INT, &new_contig);
    MPI_Type_commit(&new_contig);
    MPI_Pack(mpi_inbuf, 2, new_contig, mpi_outbuf, 20*sizeof(int), &position, MPI_COMM_WORLD);

    int res = compare_buffers(20*sizeof(int), &mpi_inbuf, &farc_inbuf, &mpi_outbuf, &farc_outbuf);
    
    free_buffers(&mpi_inbuf, &farc_inbuf, &mpi_outbuf, &farc_outbuf);
    test_result(res);

    LPK_Finalize();
    MPI_Finalize();

    return 0;

}

