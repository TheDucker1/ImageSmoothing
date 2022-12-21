/*
An implementation of the paper
"Erasing Appearance Preservation in Optimization-based Smoothing"
(European Conference on Computer Vision (ECCV) 2020)

Please check their github page for more detail (https://github.com/lllyasviel/AppearanceEraser)

Please refer below to see how the wrapper works
*/

int eap(unsigned char * image, 
    int height, int width,
    float u,
    int (*energy_smoothing_func)(unsigned char*, int, int, unsigned char* ));

/*
int L0smoothing(unsigned char * image, 
    int height, int width, int channel,
    unsigned char * mask, float kappa, float lambda);

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

    eap(image, 494, 475, 0.1, L0Smoothing_wrapper);

    //display image
    for (size_t i = 0; i < 494 * 475 * 3; i++) {
        printf("%u ", image[i]);
    }

    return 0;
}
*/