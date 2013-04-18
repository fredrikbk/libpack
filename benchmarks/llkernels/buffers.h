// Copyright 2013 Timo Schneider and Fredrik Berg Kjolstad
//
// This file is part of the libpack packing library.
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT in the top level directory for details.

#include <cstdio>
#include <cstdlib>

void init_buffers(size_t size, char* in, char* old_in, char* out, char* old_out) {
    for (size_t i=0; i<size; i++) {
        in[i] = i+1;
        old_in[i] = i+1;
        out[i] = 0;
        old_out[i] = 0;
    }

}

int compare_buffers(size_t size, char* in, char* old_in, char* out, char* old_out) {
    for (size_t i=0; i<size; i++) {
        if (in[i] != old_in[i]) return -1;
        if (out[i] != old_out[i]) return -1;
    }   

    return 0; 
}

int inspect_buffers(size_t size, char* in, char* old_in, char* out, char* old_out) {
    for (int i=0; i<size; i++) {
        printf("in[%i] = %i old_in[%i] = %i out[%i] = %i old_out[%i] = %i\n", i, in[i], i, old_in[i], i, out[i], i, old_out[i]);
    }   

    return 0; 
}

void test_result(int res) {
    if (res == 0) printf(" ... ok\n");
    else printf(" ... not ok\n");
}
