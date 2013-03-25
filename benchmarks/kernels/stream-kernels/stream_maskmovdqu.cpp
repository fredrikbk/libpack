#define NUM    8 * 1024 * 1024 
#define TYPE   double
#define ALIGN  1
#include "../stage.h"
#include <emmintrin.h>

void copy(void *out, void *in, int num, long size) {
    double *__restrict__ dout = (double*) out;
    double *__restrict__ din  = (double*) in;
    __m128d xmm1;
    __m128i mask_all=_mm_set_epi32(-1,-1,-1,-1);
    for (int i=0;i<num;i+=2) {
        xmm1 = _mm_loadu_pd(&din[i]);
        _mm_maskmoveu_si128((__m128i)xmm1, mask_all, (char*)(&dout[i]));
    }
}

