// Copyright 2013 Timo Schneider and Fredrik Berg Kjolstad
//
// This file is part of the libpack packing library.
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT in the top level directory for details.

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
