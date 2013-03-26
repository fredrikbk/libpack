#define TYPE   double
#define ALIGN  1
#include "../stage.h"
#include <cstring>

void copy(void * out, void *in, int num, long size) {
    memcpy(out, in, size);
}
