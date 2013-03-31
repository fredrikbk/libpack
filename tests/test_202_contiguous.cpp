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

    test_start("pack(2, contiguous[[int], count=0])");
    init_buffers(20*sizeof(int), &mpi_inbuf, &farc_inbuf, &mpi_outbuf, &farc_outbuf);

    farc::DDT_Init();
    farc::Datatype* t1 = new farc::PrimitiveDatatype(farc::PrimitiveDatatype::INT);
    farc::Datatype* t2 = new farc::ContiguousDatatype(t1, 0);
    farc::DDT_Commit(t2);
    farc::DDT_Pack(farc_inbuf, farc_outbuf, t2, 2);

    int position = 0;
    MPI_Datatype new_contig;
    MPI_Type_contiguous(0, MPI_INT, &new_contig);
    MPI_Type_commit(&new_contig);
    MPI_Pack(mpi_inbuf, 2, new_contig, mpi_outbuf, 20*sizeof(int), &position, MPI_COMM_WORLD);

    int res = compare_buffers(20*sizeof(int), &mpi_inbuf, &farc_inbuf, &mpi_outbuf, &farc_outbuf);
    
    free_buffers(&mpi_inbuf, &farc_inbuf, &mpi_outbuf, &farc_outbuf);
    test_result(res);

    MPI_Finalize();

    return 0;

}

