#define TYPE   double
#define ALIGN  1
#include "../stage.h"
#include <emmintrin.h>

void copy(void *out, void *in, int num, long size) {
    double *__restrict__ dout = (double*) out;
    double *__restrict__ din  = (double*) in;
    __m128d xmm1;
    for (int i=0;i<num;i+=2) {
        xmm1 = _mm_loadu_pd(&din[i]);
        _mm_storeu_pd(&dout[i], xmm1);
    }
}

