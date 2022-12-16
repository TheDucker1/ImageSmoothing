/*
An implementation of the paper
"Erasing Appearance Preservation in Optimization-based Smoothing"
(European Conference on Computer Vision (ECCV) 2020)

Please check their github page for more detail (https://github.com/lllyasviel/AppearanceEraser)

Please refer to main.c to see how the wrapper works
*/

int eap(unsigned char * image, 
    int height, int width,
    double u,
    int (*func)(unsigned char*, int, int, unsigned char* ));