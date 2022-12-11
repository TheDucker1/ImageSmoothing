#include"eap.h"
#include"L0smoothing.h"

#include<stdio.h>
#include<stdlib.h>

int L0Smoothing_wrapper(unsigned char * image, 
    int height, int width,
    unsigned char *mask) {
    return L0smoothing(image, height, width, 3, mask, 2.0, 2e-2);
}

int main() {
    //hard coded test image
    unsigned char image[494 * 475 * 3];
    for (size_t i = 0; i < 494 * 475 * 3; i++) {
        scanf("%hhu", &image[i]);
    }

    eap(image, 494, 475, 0.1, 5, L0Smoothing_wrapper);

    for (size_t i = 0; i < 494 * 475 * 3; i++) {
        printf("%u ", image[i]);
    }

    return 0;
}