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

    test_start("pack(2, [MPI_INT])");
    init_buffers(20*sizeof(int), &mpi_inbuf, &farc_inbuf, &mpi_outbuf, &farc_outbuf);

    FARC_DDT_Init();
    FARC_Datatype* t1 = new FARC_PrimitiveDatatype(FARC_PrimitiveDatatype::INT);
    FARC_DDT_Commit(t1);
    FARC_DDT_Pack(farc_inbuf, farc_outbuf, t1, 2);

    int position = 0;
    MPI_Pack(mpi_inbuf, 2, MPI_INT, mpi_outbuf, 20*sizeof(int), &position, MPI_COMM_WORLD);

    int res = compare_buffers(20*sizeof(int), &mpi_inbuf, &farc_inbuf, &mpi_outbuf, &farc_outbuf);
    free_buffers(&mpi_inbuf, &farc_inbuf, &mpi_outbuf, &farc_outbuf);
    test_result(res);

    MPI_Finalize();

    return 0;

}

