#include <string>
#include <mpi.h>

#include "test.hpp"
#include "../ddt_jit.hpp"

int main(int argc, char** argv) {

    char* mpi_inbuf;
    char* mpi_outbuf;
    char* farc_inbuf;
    char* farc_outbuf;

    struct Partstruct { 
        int foo;
        int moo[6];
        char bar[7];
    };

    struct Partstruct* particle;
    particle = (struct Partstruct*) mpi_inbuf;

    MPI_Init(&argc, &argv);

    test_start("unpack(2, struct[{1*MPI_INT, 6*MPI_INT, 8*MPI_CHAR}])");
    init_buffers(20*sizeof(int), &mpi_inbuf, &farc_inbuf, &mpi_outbuf, &farc_outbuf);

    MPI_Datatype mpitype; 
    MPI_Datatype type[3] = {MPI_INT, MPI_INT, MPI_CHAR};
    int blocklen[3] = {1, 6, 8};
    MPI_Aint disp[3];
    int base;

    // calc displacements
    MPI_Address(&(particle[0].foo), disp+0); 
    MPI_Address(&(particle[0].moo), disp+1); 
    MPI_Address(&(particle[0].bar), disp+2); 
    base = disp[0]; 
    for (int i=0; i <3; i++) disp[i] -= base;

    MPI_Type_struct(3, blocklen, disp, type, &mpitype);
    MPI_Type_commit(&mpitype);

    FARC_DDT_Init();
    FARC_Datatype* types_f[3] = {new FARC_PrimitiveDatatype(FARC_PrimitiveDatatype::INT), new FARC_PrimitiveDatatype(FARC_PrimitiveDatatype::INT), new FARC_PrimitiveDatatype(FARC_PrimitiveDatatype::CHAR)};
    FARC_Datatype* t1 = new FARC_StructDatatype(3, blocklen, disp, types_f);
    FARC_DDT_Commit(t1);
    FARC_DDT_Unpack(farc_inbuf, farc_outbuf, t1, 2);

    int position = 0;
    MPI_Unpack(mpi_inbuf, 20*sizeof(int), &position, mpi_outbuf, 2, mpitype, MPI_COMM_WORLD);

    int res = compare_buffers(20*sizeof(int), &mpi_inbuf, &farc_inbuf, &mpi_outbuf, &farc_outbuf);
    free_buffers(&mpi_inbuf, &farc_inbuf, &mpi_outbuf, &farc_outbuf);
    test_result(res);

    MPI_Finalize();

    return 0;

}

