// Copyright 2013 Timo Schneider and Fredrik Berg Kjolstad
//
// This file is part of the libpack packing library.
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT in the top level directory for details.

#include "hrtimer.h"
#include "buffers.h"

extern "C" {
    void packer(char* in, int count, char* outbuf);
    void packer_old(char* in, int count, char* outbuf);
}

unsigned long long g_timerfreq;
HRT_TIMESTAMP_T start, stop;
uint64_t packer_time, packer_old_time;

const size_t bufsize =  10 * 1024 * 1024;
char in[bufsize];
char out[bufsize] __attribute__((__aligned__(16)));
uint64_t packer_times[10];

char old_in[bufsize];
char old_out[bufsize];
uint64_t packer_old_times[10];

const int count = 2;

int main() {
    init_buffers(bufsize, in, old_in, out, old_out);
    HRT_INIT(1, g_timerfreq);
    
    // Warmup
    packer(in, count, out);
    packer(in, count, out);
    packer(in, count, out);
    packer(in, count, out);
    packer(in, count, out);
    packer(in, count, out);
    packer(in, count, out);
    packer(in, count, out);
    packer(in, count, out);
    packer(in, count, out);

    HRT_GET_TIMESTAMP(start);
    packer(in, count, out);
    HRT_GET_TIMESTAMP(stop);
    HRT_GET_ELAPSED_TICKS(start, stop, &packer_times[0]);
    HRT_GET_TIMESTAMP(start);
    packer(in, count, out);
    HRT_GET_TIMESTAMP(stop);
    HRT_GET_ELAPSED_TICKS(start, stop, &packer_times[1]);
    HRT_GET_TIMESTAMP(start);
    packer(in, count, out);
    HRT_GET_TIMESTAMP(stop);
    HRT_GET_ELAPSED_TICKS(start, stop, &packer_times[2]);
    HRT_GET_TIMESTAMP(start);
    packer(in, count, out);
    HRT_GET_TIMESTAMP(stop);
    HRT_GET_ELAPSED_TICKS(start, stop, &packer_times[3]);
    HRT_GET_TIMESTAMP(start);
    packer(in, count, out);
    HRT_GET_TIMESTAMP(stop);
    HRT_GET_ELAPSED_TICKS(start, stop, &packer_times[4]);
    HRT_GET_TIMESTAMP(start);
    packer(in, count, out);
    HRT_GET_TIMESTAMP(stop);
    HRT_GET_ELAPSED_TICKS(start, stop, &packer_times[5]);
    HRT_GET_TIMESTAMP(start);
    packer(in, count, out);
    HRT_GET_TIMESTAMP(stop);
    HRT_GET_ELAPSED_TICKS(start, stop, &packer_times[6]);
    HRT_GET_TIMESTAMP(start);
    packer(in, count, out);
    HRT_GET_TIMESTAMP(stop);
    HRT_GET_ELAPSED_TICKS(start, stop, &packer_times[7]);
    HRT_GET_TIMESTAMP(start);
    packer(in, count, out);
    HRT_GET_TIMESTAMP(stop);
    HRT_GET_ELAPSED_TICKS(start, stop, &packer_times[8]);
    HRT_GET_TIMESTAMP(start);
    packer(in, count, out);
    HRT_GET_TIMESTAMP(stop);
    HRT_GET_ELAPSED_TICKS(start, stop, &packer_times[9]);

    // Warmup
    packer_old(old_in, count, old_out);
    packer_old(old_in, count, old_out);
    packer_old(old_in, count, old_out);
    packer_old(old_in, count, old_out);
    packer_old(old_in, count, old_out);
    packer_old(old_in, count, old_out);
    packer_old(old_in, count, old_out);
    packer_old(old_in, count, old_out);
    packer_old(old_in, count, old_out);
    packer_old(old_in, count, old_out);

    HRT_GET_TIMESTAMP(start);
    packer_old(old_in, count, old_out);
    HRT_GET_TIMESTAMP(stop);
    HRT_GET_ELAPSED_TICKS(start, stop, &packer_old_times[0]);
    HRT_GET_TIMESTAMP(start);
    packer_old(old_in, count, old_out);
    HRT_GET_TIMESTAMP(stop);
    HRT_GET_ELAPSED_TICKS(start, stop, &packer_old_times[1]);
    HRT_GET_TIMESTAMP(start);
    packer_old(old_in, count, old_out);
    HRT_GET_TIMESTAMP(stop);
    HRT_GET_ELAPSED_TICKS(start, stop, &packer_old_times[2]);
    HRT_GET_TIMESTAMP(start);
    packer_old(old_in, count, old_out);
    HRT_GET_TIMESTAMP(stop);
    HRT_GET_ELAPSED_TICKS(start, stop, &packer_old_times[3]);
    HRT_GET_TIMESTAMP(start);
    packer_old(old_in, count, old_out);
    HRT_GET_TIMESTAMP(stop);
    HRT_GET_ELAPSED_TICKS(start, stop, &packer_old_times[4]);
    HRT_GET_TIMESTAMP(start);
    packer_old(old_in, count, old_out);
    HRT_GET_TIMESTAMP(stop);
    HRT_GET_ELAPSED_TICKS(start, stop, &packer_old_times[5]);
    HRT_GET_TIMESTAMP(start);
    packer_old(old_in, count, old_out);
    HRT_GET_TIMESTAMP(stop);
    HRT_GET_ELAPSED_TICKS(start, stop, &packer_old_times[6]);
    HRT_GET_TIMESTAMP(start);
    packer_old(old_in, count, old_out);
    HRT_GET_TIMESTAMP(stop);
    HRT_GET_ELAPSED_TICKS(start, stop, &packer_old_times[7]);
    HRT_GET_TIMESTAMP(start);
    packer_old(old_in, count, old_out);
    HRT_GET_TIMESTAMP(stop);
    HRT_GET_ELAPSED_TICKS(start, stop, &packer_old_times[8]);
    HRT_GET_TIMESTAMP(start);
    packer_old(old_in, count, old_out);
    HRT_GET_TIMESTAMP(stop);
    HRT_GET_ELAPSED_TICKS(start, stop, &packer_old_times[9]);


    packer_time = packer_times[0];
    packer_old_time = packer_old_times[0];
    for (int i=1;i<10;i++) {
        packer_time = (packer_times[i] < packer_time) ? packer_times[i] : packer_time;
        packer_old_time = (packer_old_times[i] < packer_old_time) ? packer_old_times[i] : packer_old_time;
    }
    printf("New: %.3lf\nOld: %.3lf\n", HRT_GET_USEC(packer_time), HRT_GET_USEC(packer_old_time));
    
    int res = compare_buffers(bufsize, in, old_in, out, old_out);
    test_result(res);
}

