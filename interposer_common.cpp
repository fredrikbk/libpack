#include <stdio.h>

extern "C" {

void interposer_init_() {
    printf("Initializing interposer\n");
}

void interposer_pack_(double *buf) {
}

}
