#define TYPE   double
#define ALIGN  1
#include "../stage.h"

void copy(void * out, void *in, int num, long size) {
    double *__restrict__ dout = (double*) out;
    double *__restrict__  din  = (double*) in;
    for (int i=0;i<num;i++) {
        dout[i] = din[i];
    }
}
