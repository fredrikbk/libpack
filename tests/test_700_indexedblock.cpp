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

    test_start("pack(2, indexed_block[blocklen=2, disp=(0,5,8), MPI_INT])");
    init_buffers(20*sizeof(int), &mpi_inbuf, &farc_inbuf, &mpi_outbuf, &farc_outbuf);

    MPI_Datatype mpitype; 
    int disp[3] = {0, 5, 8};
    int blocklen = 2;

    MPI_Type_create_indexed_block(3, blocklen, disp, MPI_INT, &mpitype);
    MPI_Type_commit(&mpitype);

    FARC_DDT_Init();
    FARC_Datatype* t1 = new FARC_PrimitiveDatatype(FARC_PrimitiveDatatype::INT);
    FARC_Datatype* t2 = new FARC_IndexedBlockDatatype(3, blocklen, disp, t1);
    FARC_DDT_Commit(t2);
    FARC_DDT_Pack(farc_inbuf, farc_outbuf, t2, 2);

    int position = 0;
    MPI_Pack(mpi_inbuf, 2, mpitype, mpi_outbuf, 20*sizeof(int), &position, MPI_COMM_WORLD);

    int res = compare_buffers(20*sizeof(int), &mpi_inbuf, &farc_inbuf, &mpi_outbuf, &farc_outbuf);

    free_buffers(&mpi_inbuf, &farc_inbuf, &mpi_outbuf, &farc_outbuf);
    test_result(res);

    MPI_Finalize();

    return 0;

}

