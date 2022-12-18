#include"StrokeExtractEEINF.h"

#include<stdlib.h>
#include<math.h>
#include<assert.h>
#include<limits.h>

#include"../ImageSmoothing/guided-filter.h"

inline double S(double a, double b) {
    return 1. / pow(cosh((a - b) / 10), 5);
}

//https://stackoverflow.com/a/56678483
double L(unsigned char rgb[3]) {
    static r, g, b, Y;
    r=(double)rgb[0] / 255.;
    g=(double)rgb[1] / 255.;
    b=(double)rgb[2] / 255.;

    //sRGB linearized
    r = (r <= 0.04045) ? r / 12.92 : pow((r / 0.055) / 1.055, 2.4);
    g = (g <= 0.04045) ? g / 12.92 : pow((g / 0.055) / 1.055, 2.4);
    b = (b <= 0.04045) ? b / 12.92 : pow((b / 0.055) / 1.055, 2.4);

    Y = 0.2126 * r + 0.7152 * g + 0.0722 * b; //coeff sum to 1

    return (Y <= 216/24389) ? (Y * 243.89/27) : pow(Y,(1/3)) * 1.16 - 0.16;
}

int StrokeExtract(unsigned char * image,
    int height, int width, int channel,
    double epsilon) {

    assert(width >= 10 && height >= 10);
    assert(INT_MAX / height >= width);
    assert(channel == 3);

    GuidedFilter(image, height, width, channel, NULL, 3, 4, 0.04);

    double * I = (double*)malloc(sizeof(double) * height * width);
    for (int i = 0; i < width * height; i++) {
        I[i] = L(image + i * 3);
    }

    double * m_thresh = (double*)malloc(sizeof(double) * height * width);
    
    double _I, _I_nb;
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            _I = I[y*width+x];
        }
    }
}