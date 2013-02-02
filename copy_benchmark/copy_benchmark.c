#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <malloc.h>

#include "hrtimer/hrtimer.h"

unsigned long long g_timerfreq;

void * sse_memcpy(void *to, const void *from, size_t len)  {
 
 unsigned long src = (unsigned long)from; 
 unsigned long dst = (unsigned long)to; 
 void *p = to; 
 int i; 
 
 /* check alignment */ 
 if ((src ^ dst) & 0xf) 
 goto unaligned; 
 
 if (src & 0xf) { 
 uint8_t chunk = 0x10 - (src & 0xf); 
 
 /* copy chunk until next 16-byte */ 
 memcpy(to, from, chunk); 
 len -= chunk; 
 to += chunk; 
 from += chunk; 
 } 
 
 /* 
 * copy in 256 Byte portions 
 */ 
 for (i = 0; i < (len & ~0xff); i += 256) { 
 asm volatile(
   "prefetchnta 128(%%rsi)\n\t" //SSE2 prefetch
   "prefetchnta 160(%%rsi)\n\t"
   "prefetchnta 192(%%rsi)\n\t"
   "prefetchnta 224(%%rsi)\n\t"
 
 "movaps 0x0(%0), %%xmm0\n\t" 
 "movaps 0x10(%0), %%xmm1\n\t" 
 "movaps 0x20(%0), %%xmm2\n\t" 
 "movaps 0x30(%0), %%xmm3\n\t" 
 "movaps 0x40(%0), %%xmm4\n\t" 
 "movaps 0x50(%0), %%xmm5\n\t" 
 "movaps 0x60(%0), %%xmm6\n\t" 
 "movaps 0x70(%0), %%xmm7\n\t" 
 "movaps 0x80(%0), %%xmm8\n\t" 
 "movaps 0x90(%0), %%xmm9\n\t" 
 "movaps 0xa0(%0), %%xmm10\n\t" 
 "movaps 0xb0(%0), %%xmm11\n\t" 
 "movaps 0xc0(%0), %%xmm12\n\t" 
 "movaps 0xd0(%0), %%xmm13\n\t" 
 "movaps 0xe0(%0), %%xmm14\n\t" 
 "movaps 0xf0(%0), %%xmm15\n\t" 
 
 "movaps %%xmm0, 0x0(%1)\n\t" 
 "movaps %%xmm1, 0x10(%1)\n\t" 
 "movaps %%xmm2, 0x20(%1)\n\t" 
 "movaps %%xmm3, 0x30(%1)\n\t" 
 "movaps %%xmm4, 0x40(%1)\n\t" 
 "movaps %%xmm5, 0x50(%1)\n\t" 
 "movaps %%xmm6, 0x60(%1)\n\t" 
 "movaps %%xmm7, 0x70(%1)\n\t" 
 "movaps %%xmm8, 0x80(%1)\n\t" 
 "movaps %%xmm9, 0x90(%1)\n\t" 
 "movaps %%xmm10, 0xa0(%1)\n\t" 
 "movaps %%xmm11, 0xb0(%1)\n\t" 
 "movaps %%xmm12, 0xc0(%1)\n\t" 
 "movaps %%xmm13, 0xd0(%1)\n\t" 
 "movaps %%xmm14, 0xe0(%1)\n\t" 
 "movaps %%xmm15, 0xf0(%1)\n\t" 
 : : "r" (from), "r" (to) : "memory"); 
 
 from += 256; 
 to += 256; 
 } 
 
 goto trailer; 
 
unaligned: 
 /* 
 * copy in 256 Byte portions unaligned 
 */ 
 for (i = 0; i < (len & ~0xff); i += 256) { 
 asm volatile( 
 "movups 0x0(%0), %%xmm0\n\t" 
 "movups 0x10(%0), %%xmm1\n\t" 
 "movups 0x20(%0), %%xmm2\n\t" 
 "movups 0x30(%0), %%xmm3\n\t" 
 "movups 0x40(%0), %%xmm4\n\t" 
 "movups 0x50(%0), %%xmm5\n\t" 
 "movups 0x60(%0), %%xmm6\n\t" 
 "movups 0x70(%0), %%xmm7\n\t" 
 "movups 0x80(%0), %%xmm8\n\t" 
 "movups 0x90(%0), %%xmm9\n\t" 
 "movups 0xa0(%0), %%xmm10\n\t" 
 "movups 0xb0(%0), %%xmm11\n\t" 
 "movups 0xc0(%0), %%xmm12\n\t" 
 "movups 0xd0(%0), %%xmm13\n\t" 
 "movups 0xe0(%0), %%xmm14\n\t" 
 "movups 0xf0(%0), %%xmm15\n\t" 
 
 "movups %%xmm0, 0x0(%1)\n\t" 
 "movups %%xmm1, 0x10(%1)\n\t" 
 "movups %%xmm2, 0x20(%1)\n\t" 
 "movups %%xmm3, 0x30(%1)\n\t" 
 "movups %%xmm4, 0x40(%1)\n\t" 
 "movups %%xmm5, 0x50(%1)\n\t" 
 "movups %%xmm6, 0x60(%1)\n\t" 
 "movups %%xmm7, 0x70(%1)\n\t" 
 "movups %%xmm8, 0x80(%1)\n\t" 
 "movups %%xmm9, 0x90(%1)\n\t" 
 "movups %%xmm10, 0xa0(%1)\n\t" 
 "movups %%xmm11, 0xb0(%1)\n\t" 
 "movups %%xmm12, 0xc0(%1)\n\t" 
 "movups %%xmm13, 0xd0(%1)\n\t" 
 "movups %%xmm14, 0xe0(%1)\n\t" 
 "movups %%xmm15, 0xf0(%1)\n\t" 
 : : "r" (from), "r" (to) : "memory"); 
 
 from += 256; 
 to += 256; 
 } 
 
trailer: 
 memcpy(to, from, len & 0xff); 
 
 return p; 
} 

double memcpy_benchmark(size_t size, char* src, char* dst) {

    HRT_TIMESTAMP_T start, stop;
    uint64_t elapsed_ticks;

    HRT_GET_TIMESTAMP(start);
    memcpy(dst, src, size);
    HRT_GET_TIMESTAMP(stop);

    HRT_GET_ELAPSED_TICKS(start, stop, &elapsed_ticks);

    return HRT_GET_USEC(elapsed_ticks);

}

double assignment_benchmark(size_t size, char* src, char* dst) {

    HRT_TIMESTAMP_T start, stop;
    uint64_t elapsed_ticks;

    HRT_GET_TIMESTAMP(start);
    for (int i=0; i<size; i++) {
        dst[i] = src[i];
    }
    HRT_GET_TIMESTAMP(stop);

    HRT_GET_ELAPSED_TICKS(start, stop, &elapsed_ticks);

    return HRT_GET_USEC(elapsed_ticks);

}

double prefetch_benchmark(size_t size, char* src, char* dst) {

    HRT_TIMESTAMP_T start, stop;
    uint64_t elapsed_ticks;

    HRT_GET_TIMESTAMP(start);

    int dist = 32;
    for (int i=0; i<size; i++) {
        dst[i] = src[i];
        __builtin_prefetch (&src[i+dist], 0, 0 /*reading, no temp. locality*/);
        __builtin_prefetch (&dst[i+dist], 1, 0 /*writing, no temp. locality*/);
    }
    HRT_GET_TIMESTAMP(stop);

    HRT_GET_ELAPSED_TICKS(start, stop, &elapsed_ticks);

    return HRT_GET_USEC(elapsed_ticks);

}

void memcpy2(void* dest, const void* src, const unsigned long size_t) {

    unsigned long size1, size2 = 0, size3;

    // can we synchronize src and dst pointers?
    if (((((unsigned long) dest) ^ ((unsigned long) src)) & 0xf) == 0) {
        // if yes, compute the part1 and part 2 - part 2 is the count of 128 byte block
        size1 = (0x10 - (((unsigned long) dest) & 0xf)) & 0xf;
        size2 = (size_t - size1) & ~0x7F;
    }

  // if the size that can be copied by 128 blocks (size2) is > 0 then use fast copy
  if (size2) {     
    if (size1) memcpy(dest, src, size1);
  
    asm("push %rsi\n\t");                      // preserve ESI to the stack
    asm("push %rdi\n\t");                      // preserve EDI to the stack
  
    asm("mov %0, %%rsi\n\t"
      : :"r"(src));     //src pointer

    asm("mov %0, %%rdi\n\t"
      :
      : "r"(dest));  //dest pointer

    asm("mov %0, %%rbx\n\t"  //ebx is our counter
      :
      : "r"(size2));

    asm("add %0, %%rsi\n\t" :: "r"(size1));     // add the copied bytes count to src pointer
    asm("add %0, %%rdi\n\t" :: "r"(size1));     // add the copied bytes count to dest pointer

    asm("shr $7, %rbx\n\t"      //divide by 128 (8 * 128bit registers)

      "loop_copy:\n\t"
        "prefetchnta 128(%rsi)\n\t" //SSE2 prefetch
        "prefetchnta 160(%rsi)\n\t"
        "prefetchnta 192(%rsi)\n\t"
        "prefetchnta 224(%rsi)\n\t"

        "movdqa 0(%rsi), %xmm0\n\t" //move data from src to registers
        "movdqa 16(%rsi), %xmm1\n\t"
        "movdqa 32(%rsi), %xmm2\n\t"
        "movdqa 48(%rsi), %xmm3\n\t"
        "movdqa 64(%rsi), %xmm4\n\t"
        "movdqa 80(%rsi), %xmm5\n\t"
        "movdqa 96(%rsi), %xmm6\n\t"
        "movdqa 112(%rsi), %xmm7\n\t"
 
        "movntdq %xmm0, 0(%rdi)\n\t" //move data from registers to dest
        "movntdq %xmm1, 16(%rdi)\n\t"
        "movntdq %xmm2, 32(%rdi)\n\t"
        "movntdq %xmm3, 48(%rdi)\n\t"
        "movntdq %xmm4, 64(%rdi)\n\t"
        "movntdq %xmm5, 80(%rdi)\n\t"
        "movntdq %xmm6, 96(%rdi)\n\t"
        "movntdq %xmm7, 112(%rdi)\n\t"

        "add $128, %rsi\n"
        "add $128, %rdi\n"
        "dec %rbx\n"

        "jnz loop_copy\n" //loop please
      "loop_copy_end:\n"
    );
 
    asm ("pop %rdi\n\t");                      // restore EDI
    asm ("pop %rsi\n\t");                      // restore ESI
  
    size3 = size_t - size2 - size1;
    if (size3) 
      memcpy(&(((char*)dest)[size1 + size2]), 
        &(((char*) src)[size1 + size2]), size3);
  }
  else
  {
      printf("memcpy used\n");
      memcpy(dest,src,size_t);
  }
}

double sse_benchmark(size_t size, char* src, char* dst) {

    HRT_TIMESTAMP_T start, stop;
    uint64_t elapsed_ticks;

    HRT_GET_TIMESTAMP(start);

    sse_memcpy(dst, src, size);    

    HRT_GET_TIMESTAMP(stop);

    HRT_GET_ELAPSED_TICKS(start, stop, &elapsed_ticks);

    return HRT_GET_USEC(elapsed_ticks);

}

void prep_benchmarks(size_t size, char** src, char** dst) {

    *src = memalign(128, size);
    *dst = memalign(128, size);

    assert(src != NULL);
    assert(dst != NULL);

    // avoid overoptimization
    for (int i=0; i<size; i++) (*src)[i] = rand() % 42;
    for (int i=0; i<size; i++) (*dst)[i] = 1;

}

void fin_benchmarks(int size, char* src, char* dst) {

    free(src);
    free(dst);    

    // avoid overoptimization
    uint64_t res = 0;
    for (int i=0; i<size; i++) res += dst[i];
    if (res == 0) printf("");

}

int main(int argc, char **argv) {
    
    char *src, *dst; 
    HRT_INIT(1, g_timerfreq);

    printf("%10s %10s %10s\n", "size", "method", "time");

    for (int size=1; size<1e8; size*=2) {

       prep_benchmarks(size, &src, &dst);

        for (int run=0; run<25; run++) {
            double time_assignment = assignment_benchmark(size, src, dst);
            printf("%10i %10s %10.3lf\n", size, "assign", time_assignment);
        }

        for (int run=0; run<25; run++) {
            double time_memcpy = memcpy_benchmark(size, src, dst);
            printf("%10i %10s %10.3lf\n", size, "memcpy", time_memcpy);
        }

        for (int run=0; run<25; run++) {
            double time_prefetch = prefetch_benchmark(size, src, dst);
            printf("%10i %10s %10.3lf\n", size, "prefetch", time_prefetch);
        }

        for (int run=0; run<25; run++) {
            double time_sse = sse_benchmark(size, src, dst);
            printf("%10i %10s %10.3lf\n", size, "sse", time_sse);
        }

        fin_benchmarks(size, src, dst);

        
    }

}
