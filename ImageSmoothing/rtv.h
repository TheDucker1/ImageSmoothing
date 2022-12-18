/*
An implementation of the paper
Structure Extraction from Texture via Relative Total Variation
(The Chinese Univeristy of Hong Kong)

Based on the implemetation of the ToS2P (https://github.com/lllyasviel)

*/

int rtv(unsigned char * image, 
    int height, int width, int channel,
    unsigned char * mask,
    double lambda, double sigma, double sharpness, int maxIter);