// Copyright 2013 Timo Schneider and Fredrik Berg Kjolstad
//
// This file is part of the libpack packing library.
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT in the top level directory for details.

#include "hrtimer.h"

void stage(char **out, char **in, int *num, long *size);
void  copy(void *out,  void  *in, int  num, long  size);

unsigned long long g_timerfreq;
HRT_TIMESTAMP_T start, stop;
uint64_t copy_time;
uint64_t copy_times[10];

int main(int argc, char *argv[]) {
    HRT_INIT(0, g_timerfreq);
    
    char *in;
    char *out;
    int  num;
    long size;
    stage(&out, &in, &num, &size);

    for (int i=0; i<num; i++) {
        in[i] = i;
        out[i] = 0;
    }
    
    // Warmup
    copy(out, in, num, size);
    copy(out, in, num, size);
    copy(out, in, num, size);
    copy(out, in, num, size);
    copy(out, in, num, size);

    HRT_GET_TIMESTAMP(start);
    copy(out, in, num, size);
    HRT_GET_TIMESTAMP(stop);
    HRT_GET_ELAPSED_TICKS(start, stop, &copy_times[0]);
    HRT_GET_TIMESTAMP(start);
    copy(out, in, num, size);
    HRT_GET_TIMESTAMP(stop);
    HRT_GET_ELAPSED_TICKS(start, stop, &copy_times[1]);
    HRT_GET_TIMESTAMP(start);
    copy(out, in, num, size);
    HRT_GET_TIMESTAMP(stop);
    HRT_GET_ELAPSED_TICKS(start, stop, &copy_times[2]);
    HRT_GET_TIMESTAMP(start);
    copy(out, in, num, size);
    HRT_GET_TIMESTAMP(stop);
    HRT_GET_ELAPSED_TICKS(start, stop, &copy_times[3]);
    HRT_GET_TIMESTAMP(start);
    copy(out, in, num, size);
    HRT_GET_TIMESTAMP(stop);
    HRT_GET_ELAPSED_TICKS(start, stop, &copy_times[4]);
    HRT_GET_TIMESTAMP(start);
    copy(out, in, num, size);
    HRT_GET_TIMESTAMP(stop);
    HRT_GET_ELAPSED_TICKS(start, stop, &copy_times[5]);
    HRT_GET_TIMESTAMP(start);
    copy(out, in, num, size);
    HRT_GET_TIMESTAMP(stop);
    HRT_GET_ELAPSED_TICKS(start, stop, &copy_times[6]);
    HRT_GET_TIMESTAMP(start);
    copy(out, in, num, size);
    HRT_GET_TIMESTAMP(stop);
    HRT_GET_ELAPSED_TICKS(start, stop, &copy_times[7]);
    HRT_GET_TIMESTAMP(start);
    copy(out, in, num, size);
    HRT_GET_TIMESTAMP(stop);
    HRT_GET_ELAPSED_TICKS(start, stop, &copy_times[8]);
    HRT_GET_TIMESTAMP(start);
    copy(out, in, num, size);
    HRT_GET_TIMESTAMP(stop);
    HRT_GET_ELAPSED_TICKS(start, stop, &copy_times[9]);

    copy_time = copy_times[0];
    for (int i=1;i<10;i++) {
        copy_time = (copy_times[i] < copy_time) ? copy_times[i] : copy_time;
    }
    double copy_time_usec = HRT_GET_USEC(copy_time);
    printf("%19s:  %8.2lf us   %6.3lf (%6.3lf)GB/s    (%d doubles, %d bytes)\n", argv[0], copy_time_usec,
            ((size) / (copy_time_usec*1024*1024*1024)*1000000),
            ((size) / (HRT_GET_USEC(copy_times[0])*1024*1024*1024)*1000000),
            num, size);

    for (int i=0; i<num; i++) {
        if (out[i] != in[i]) {
            printf("ERROR (%d): %lf != %lf\n", i, out[i], in[i]);
            return 1;
        }
    }
}

