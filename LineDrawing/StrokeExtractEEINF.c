#include"StrokeExtractEEINF.h"

#include<stdlib.h>
#include<math.h>
#include<assert.h>
#include<limits.h>

//#include<omp.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define util_max(a, b) (((a) > (b)) ? (a) : (b))
#define util_min(a, b) (((a) < (b)) ? (a) : (b))
#define util_square(a) ((a)*(a))
#define util_inv(a) (1. / (a))  //imply double

inline double S(double a, double b) {
    return 1. / pow(cosh((a - b) / 10), 5);
}

//https://stackoverflow.com/a/56678483
double L(unsigned char rgb[3]) {
    static double r, g, b, Y;
    r=(double)rgb[0] / 255.;
    g=(double)rgb[1] / 255.;
    b=(double)rgb[2] / 255.;

    //sRGB linearized
    r = (r <= 0.04045) ? r / 12.92 : pow((r + 0.055) / 1.055, 2.4);
    g = (g <= 0.04045) ? g / 12.92 : pow((g + 0.055) / 1.055, 2.4);
    b = (b <= 0.04045) ? b / 12.92 : pow((b + 0.055) / 1.055, 2.4);

    Y = 0.2126 * r + 0.7152 * g + 0.0722 * b; //coeff sum to 1
    //Y = (Y <= 0.008856) ? Y * (903.3) : pow(Y,(1/3)) * 116 - 16;
    //Y /= 100;

    return Y;
    //range [0, 1]
}

int StrokeExtract(unsigned char * image,
    int height, int width, int channel,
    int radius, double epsilon) {

    assert(width >= 10 && height >= 10);
    assert(INT_MAX / height >= width);
    assert(channel == 3);

    //GuidedFilter(image, height, width, channel, NULL, 3, 4, 0.04);
    //assuming the image was smoothed beforehand

    double * I = (double*)malloc(sizeof(double) * height * width);
    for (int i = 0; i < width * height; i++) {
        I[i] = L(image + i * 3);
    }

    //threshold mass of similarity
    double * m_thresh = (double*)calloc(height * width, sizeof(double));

    const int radius_square = util_square(radius);
    double _I, _I_nb;
    int y, x, dy, dx, ny, nx;
    
    //#pragma omp parallel for
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            _I = I[y*width+x];
            if (_I > epsilon) {
                m_thresh[y*width+x] = 1;
                continue;
            }
            for (dy = -radius; dy <= radius; dy++) {
                ny = y + dy;
                if (ny < 0) ny = 0;
                else if (ny >= height) ny = height - 1;

                for (dx = -radius; dx <= radius; dx++) {
                    nx = x + dx;
                    if (dx * dx + dy * dy > radius_square) {
                        continue;
                    } 
                    if (nx < 0) nx = 0;
                    else if (nx >= width) nx = width - 1;

                    _I_nb = I[ny*width + nx];
                    m_thresh[y*width + x] += S(_I, _I_nb);
                    
                }
            }
        }
    }

    //LoG
    //https://stackoverflow.com/a/53665075
    double * LoG = (double*)calloc(height * width, sizeof(double));
    double * tmp1 = (double*)malloc(sizeof(double) * height * width);

    const double sigma = 1.1;
    const int cutoff = ceil(sigma * 3);
    const double sigma_pow2 = util_square(sigma);
    const double sigma_pow4 = util_square(sigma_pow2);
    int ksize = 2 * cutoff + 1;
    int r = ksize / 2;
    double * kernel = (double*)malloc(sizeof(double) * ksize), s = 0;
    double * Lkernel = (double*)malloc(sizeof(double) * ksize);
    for (int i = 0; i <= r; i++) {
        kernel[r-i] = kernel[r+i] = util_inv(exp((double)(i*i) / (2 * (sigma*sigma)))) * util_inv(sigma * sqrt(2 * M_PI));
        if (i == 0) s += kernel[r-i];
        else s += 2 * kernel[r-i];
    }
    for (int i = 0; i < ksize; i++) {
        kernel[i] /= s;
    }
    for (int i = -cutoff; i <= cutoff; i++) {
        Lkernel[i+cutoff] = kernel[i+cutoff] * (util_square(i) - sigma_pow2) / sigma_pow4;
    }

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            tmp1[y*width + x] = 0;
            for (int dx = -cutoff; dx <= cutoff; dx++) {
                int nx = x + dx;
                if (nx < 0) nx = 0;
                else if (nx >= width) nx = width - 1;

                tmp1[y*width + x] += kernel[dx+cutoff] * I[y * width + nx];
            }
        }
    }
    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            for (int dy = -cutoff; dy <= cutoff; dy++) {
                int ny = y + dy;
                if (ny < 0) ny = 0;
                else if (ny >= height) ny = height - 1;

                LoG[y*width + x] += Lkernel[dy+cutoff] * tmp1[ny*width + x];
            }
        }
    }

    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            tmp1[y*width + x] = 0;
            for (int dy = -cutoff; dy <= cutoff; dy++) {
                int ny = y + dy;
                if (ny < 0) ny = 0;
                else if (ny >= height) ny = height - 1;

                tmp1[y*width + x] += kernel[dy+cutoff] * I[ny*width + x];
            }
        }
    }
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            for (int dx = -cutoff; dx <= cutoff; dx++) {
                int nx = x + dx;
                if (nx < 0) nx = 0;
                else if (nx >= width) nx = width - 1;

                LoG[y*width + x] += Lkernel[dx+cutoff] * tmp1[y * width + nx];
            }
            _I = LoG[y*width + x];
           LoG[y*width + x] = (_I > 0) ? _I : 0;
        }
    }

    free(kernel); free(Lkernel);
    free(tmp1);

    double _m_edge;
    const double g = M_PI * util_square(r) / 2.;
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            _m_edge = LoG[y*width + x] * m_thresh[y*width + x];
            I[y*width + x] = (_m_edge <= g) ? 1 - _m_edge / g : 0;
        }
    }

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            _I = I[y*width+x];
            image[3 * (y*width+x)] = image[3 * (y*width+x) + 1] = image[3 * (y*width+x) + 2] = (_I < 0) ? 0 : ((_I > 1) ? 255 : (unsigned char)(_I * 255));
        }
    }


    free(LoG);
    free(m_thresh);
    free(I);

    return 0;
}