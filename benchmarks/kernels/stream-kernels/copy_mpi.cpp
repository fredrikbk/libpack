#define TYPE   double
#define ALIGN  1
#include <mpi.h>

MPI_Datatype ddt;

TYPE  inbuf[NUM] __attribute__((aligned(ALIGN)));
TYPE outbuf[NUM] __attribute__((aligned(ALIGN)));
void stage(char **out, char **in, int *num, long *size) {
    char **argv;
    int argc = 0;
    MPI_Init(&argc, &argv);
    ddt = MPI_DOUBLE;

    *in   = (char*)inbuf;
    *out  = (char*)outbuf;
    *num  = NUM;
    *size = NUM*sizeof(TYPE);
}

void copy(void * out, void *in, int num, long size) {
    int position = 0;
    MPI_Pack(in, num, ddt, out, size, &position, MPI_COMM_WORLD);
}
