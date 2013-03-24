#include <string>
#include <mpi.h>

#include "../ddt_jit.hpp"
#include "test.hpp"

#define count    1
#define blocklen 7
#define stride   8

int main(int argc, char** argv) {

    int position = 0;
    char* mpi_inbuf;
    char* mpi_outbuf;
    char* farc_inbuf;
    char* farc_outbuf;

    MPI_Init(&argc, &argv);

    test_start("vector_bug");

    farc::DDT_Init();
    farc::Datatype* t1 = new farc::PrimitiveDatatype(farc::PrimitiveDatatype::DOUBLE);
    farc::Datatype* t2 = new farc::VectorDatatype(t1, count, blocklen, stride);
    farc::DDT_Commit(t2);

    MPI_Datatype newtype;
    MPI_Type_vector(count, blocklen, stride, MPI_DOUBLE, &newtype);
    MPI_Type_commit(&newtype);

    int extent = t2->getExtent();
    init_buffers(extent, &mpi_inbuf, &farc_inbuf, &mpi_outbuf, &farc_outbuf);

    farc::DDT_Pack(farc_inbuf, farc_outbuf, t2, 1);
    MPI_Pack(mpi_inbuf, 1, newtype, mpi_outbuf, extent, &position, MPI_COMM_WORLD);

    int res = compare_buffers(extent, &mpi_inbuf, &farc_inbuf, &mpi_outbuf, &farc_outbuf);
    free_buffers(&mpi_inbuf, &farc_inbuf, &mpi_outbuf, &farc_outbuf);
    test_result(res);

    MPI_Finalize();

    return 0;

}

