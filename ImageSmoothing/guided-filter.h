/*
Implementation of the paper
"Guided Image Filtering"
ECCV 2010

Based on the original authors' implementation in Matlab
*/

int GuidedFilter(unsigned char * image,
    int height, int width, int channel,
    unsigned char * guidance,
    int guidance_channel,
    int r, float epsilon);