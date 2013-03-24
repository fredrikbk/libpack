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

    test_start("free subtype before usage of constructed type");
    init_buffers(20*sizeof(int), &mpi_inbuf, &farc_inbuf, &mpi_outbuf, &farc_outbuf);

    farc::DDT_Init();
    farc::Datatype* t1 = new farc::PrimitiveDatatype(farc::PrimitiveDatatype::INT);
    farc::Datatype* t2 = new farc::ContiguousDatatype(t1, 4);
    farc::Datatype* t3 = new farc::ContiguousDatatype(t2, 2);
    farc::DDT_Free(t2);
    farc::DDT_Commit(t3);
    farc::DDT_Pack(farc_inbuf, farc_outbuf, t3, 2);

    int position = 0;
    MPI_Datatype new_contig;
    MPI_Datatype newer_contig;
    MPI_Type_contiguous(4, MPI_INT, &new_contig);
    MPI_Type_contiguous(2, new_contig, &newer_contig);
    MPI_Type_free(&new_contig);
    MPI_Type_commit(&newer_contig);
    MPI_Pack(mpi_inbuf, 2, newer_contig, mpi_outbuf, 20*sizeof(int), &position, MPI_COMM_WORLD);

    int res = compare_buffers(20*sizeof(int), &mpi_inbuf, &farc_inbuf, &mpi_outbuf, &farc_outbuf);
    
    free_buffers(&mpi_inbuf, &farc_inbuf, &mpi_outbuf, &farc_outbuf);
    test_result(res);

    MPI_Finalize();

    return 0;

}

