TYPE  inbuf[NUM] __attribute__((aligned(ALIGN)));
TYPE outbuf[NUM] __attribute__((aligned(ALIGN)));
void stage(char **out, char **in, int *num, long *size) {
    *in   = (char*)inbuf;
    *out  = (char*)outbuf;
    *num  = NUM;
    *size = NUM*sizeof(TYPE);
}
