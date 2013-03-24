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
        char bar[8];
    };

    struct Partstruct* particle;
    particle = (struct Partstruct*) mpi_inbuf;

    MPI_Init(&argc, &argv);


    test_start("pack(2, struct[{1*MPI_INT, 6*MPI_INT, 8*MPI_CHAR}])");
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

    farc::DDT_Init();
    farc::Datatype* types_f[3] = {new farc::PrimitiveDatatype(farc::PrimitiveDatatype::INT), new farc::PrimitiveDatatype(farc::PrimitiveDatatype::INT), new farc::PrimitiveDatatype(farc::PrimitiveDatatype::CHAR)};
    farc::Datatype* t1 = new farc::StructDatatype(3, blocklen, disp, types_f);
    farc::DDT_Commit(t1);
    farc::DDT_Pack(farc_inbuf, farc_outbuf, t1, 2);

    int position = 0;
    MPI_Pack(mpi_inbuf, 2, mpitype, mpi_outbuf, 20*sizeof(int), &position, MPI_COMM_WORLD);

    int res = compare_buffers(20*sizeof(int), &mpi_inbuf, &farc_inbuf, &mpi_outbuf, &farc_outbuf);
    free_buffers(&mpi_inbuf, &farc_inbuf, &mpi_outbuf, &farc_outbuf);
    test_result(res);

    MPI_Finalize();

    return 0;

}

