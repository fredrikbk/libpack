#define NUM    8 * 1024 * 1024 
#define TYPE   double
#define ALIGN  1
#include "../stage.h"
#include <immintrin.h>

void copy(void *out, void *in, int num, long size) {
    double *__restrict__ dout = (double*) out;
    double *__restrict__ din  = (double*) in;
    __m256d ymm1;
    for (int i=0;i<num;i+=4) {
        ymm1 = _mm256_loadu_pd(&din[i]);
        _mm256_storeu_pd(&dout[i], ymm1);
    }
}

