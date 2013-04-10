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

    test_start("pack(2, hvector[[int], count=2, blklen=3, stride=5*sizeof(int)])");
    init_buffers(20*sizeof(int), &mpi_inbuf, &farc_inbuf, &mpi_outbuf, &farc_outbuf);

    farc::DDT_Init();
    farc::Datatype* t1 = new farc::PrimitiveDatatype(farc::PrimitiveDatatype::INT);
    farc::Datatype* t2 = new farc::HVectorDatatype(2, 3, 5*sizeof(int), t1);
    farc::DDT_Commit(t2);
    farc::DDT_Pack(farc_inbuf, farc_outbuf, t2, 2);

    MPI_Datatype newtype;
    MPI_Type_hvector(2, 3, 5*sizeof(int), MPI_INT, &newtype);
    MPI_Type_commit(&newtype);
    int position = 0;
    MPI_Pack(mpi_inbuf, 2, newtype, mpi_outbuf, 20*sizeof(int), &position, MPI_COMM_WORLD);

    int res = compare_buffers(20*sizeof(int), &mpi_inbuf, &farc_inbuf, &mpi_outbuf, &farc_outbuf);
    free_buffers(&mpi_inbuf, &farc_inbuf, &mpi_outbuf, &farc_outbuf);
    test_result(res);

    MPI_Finalize();

    return 0;

}

