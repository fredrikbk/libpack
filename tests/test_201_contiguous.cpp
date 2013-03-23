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

    test_start("unpack(2, contiguous[[int], count=8])");
    init_buffers(20*sizeof(int), &mpi_inbuf, &farc_inbuf, &mpi_outbuf, &farc_outbuf);

    FARC_DDT_Init();
    FARC_Datatype* t1 = new FARC_PrimitiveDatatype(FARC_PrimitiveDatatype::INT);
    FARC_Datatype* t2 = new FARC_ContiguousDatatype(t1, 8);
    FARC_DDT_Commit(t2);
    FARC_DDT_Unpack(farc_inbuf, farc_outbuf, t2, 2);

    int position = 0;
    MPI_Datatype new_contig;
    MPI_Type_contiguous(8, MPI_INT, &new_contig);
    MPI_Type_commit(&new_contig);
    MPI_Unpack(mpi_inbuf, 20*sizeof(int), &position, mpi_outbuf, 2, new_contig, MPI_COMM_WORLD);

    int res = compare_buffers(20*sizeof(int), &mpi_inbuf, &farc_inbuf, &mpi_outbuf, &farc_outbuf);
    
    free_buffers(&mpi_inbuf, &farc_inbuf, &mpi_outbuf, &farc_outbuf);
    test_result(res);

    MPI_Finalize();

    return 0;

}

