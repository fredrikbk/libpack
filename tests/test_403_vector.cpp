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

    test_start("unpack(2, vector[[int], count=4, blklen=3, stride=5]) [unrolling]");
    init_buffers(40*sizeof(int), &mpi_inbuf, &farc_inbuf, &mpi_outbuf, &farc_outbuf);

    FARC_DDT_Init();
    FARC_Datatype* t1 = new FARC_PrimitiveDatatype(FARC_PrimitiveDatatype::INT);
    FARC_Datatype* t2 = new FARC_VectorDatatype(t1, 4, 3, 5);
    FARC_DDT_Commit(t2);
    FARC_DDT_Unpack(farc_inbuf, farc_outbuf, t2, 2);

    MPI_Datatype newtype;
    MPI_Type_vector(4, 3, 5, MPI_INT, &newtype);
    MPI_Type_commit(&newtype);
    int position = 0;
    MPI_Unpack(mpi_inbuf, 40*sizeof(int), &position, mpi_outbuf, 2, newtype, MPI_COMM_WORLD);

    int res = compare_buffers(40*sizeof(int), &mpi_inbuf, &farc_inbuf, &mpi_outbuf, &farc_outbuf);
    free_buffers(&mpi_inbuf, &farc_inbuf, &mpi_outbuf, &farc_outbuf);
    test_result(res);

    MPI_Finalize();

    return 0;

}

