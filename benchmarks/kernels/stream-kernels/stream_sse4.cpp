// Copyright 2013 Timo Schneider and Fredrik Berg Kjolstad
//
// This file is part of the libpack packing library.
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT in the top level directory for details.

#define TYPE   double
#define ALIGN  16
#include "../stage.h"
#include <smmintrin.h>

void copy(void *out, void *in, int num, long size) {
    double *__restrict__ dout = (double*) out;
    double *__restrict__ din  = (double*) in;
    __m128d xmm1;
    for (int i=0;i<num;i+=2) {
        xmm1 = (__m128d)_mm_stream_load_si128((__m128i*)(&din[i]));
        _mm_stream_pd(&dout[i], xmm1);
    }
}

