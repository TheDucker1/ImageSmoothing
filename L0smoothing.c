#include"L0smoothing.h"

#include<stdlib.h>
#include<math.h>
#include<assert.h>

#include"fftw3.h"

int L0smoothing_L(unsigned char * image, 
    int height, int width,
    unsigned char * mask, double kappa, double lambda) {

    assert(width > 0 && height > 0);

    double *Im = (double*)malloc(sizeof(double) * height * width * 3);
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            Im[(y*width+x)] = (double)(image[(y*width+x)]) / 255.;
        }
    }
    
    double *Ma = (double*)malloc(sizeof(double) * height * width);
    for (size_t i = 0; i < (size_t)(width * height); i++) {
        Ma[i] = (mask == NULL) ? 0. : ((double)mask[i] / 255.);
    }
    
    fftw_complex *fx_in, *fx_out;
    fftw_complex *fy_in, *fy_out;
    fftw_plan pfx, pfy;

    fx_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * width * height);
    fx_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * width * height);
    fy_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * width * height);
    fy_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * width * height);

    pfx = fftw_plan_dft_2d(height, width, fx_in, fx_out, FFTW_FORWARD, FFTW_ESTIMATE);
    pfy = fftw_plan_dft_2d(height, width, fy_in, fy_out, FFTW_FORWARD, FFTW_ESTIMATE);

    fx_in[0][0] = -1.;
    fx_in[0][1] = 0;
    fx_in[width-1][0] = 1.;
    fx_in[width-1][1] = 0;
    fy_in[0][0] = -1.;
    fy_in[0][1] = 0;
    fy_in[(height-1)*width][0] = 1.;
    fy_in[(height-1)*width][1] = 0;

    fftw_execute(pfx);
    fftw_execute(pfy);

    double *Denormin2 = (double*)malloc(sizeof(double) * width * height);
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            Denormin2[y*width + x] = pow(fx_out[y*width+x][0],2) + pow(fx_out[y*width+x][1],2)
                + pow(fy_out[y*width+x][0],2) + pow(fy_out[y*width+x][1],2);
        }
    }

    fftw_destroy_plan(pfx);
    fftw_destroy_plan(pfy);
    fftw_free(fx_in); fftw_free(fx_out);
    fftw_free(fy_in); fftw_free(fy_out);

    fftw_complex *Normin1, *Normin1In;
    fftw_plan p;
    Normin1In = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * width * height);
    Normin1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * width * height);
    p = fftw_plan_dft_2d(height, width, Normin1In, Normin1, FFTW_FORWARD, FFTW_ESTIMATE);
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            Normin1In[y*width+x][0] = Im[(y*width+x)];
            Normin1In[y*width+x][1] = 0;
        }
    }
    fftw_execute(p);
        
    fftw_destroy_plan(p);
    fftw_free(Normin1In);

    double *h = (double*)malloc(sizeof(double) * height * width);
    double *v = (double*)malloc(sizeof(double) * height * width);
    
    const double betamax = 1e5;
    double beta = 2 * lambda;

    double *Denormin = (double*)malloc(sizeof(double) * height * width);
    fftw_complex *FS, *FSIn, *Normin2, *Normin2In;
    fftw_plan FSp;
    FSIn = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * width * height);
    Normin2In = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * width * height);
    FS = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * width * height);
    Normin2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * width * height);
    
    while (beta < betamax) {
        double _h, _v;
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                Denormin[y*width+x] = 1. + beta * Denormin2[y*width+x];

                if (Ma[(y*width+x)] > 0.5) {
                    h[(y*width+x)] = 0;
                    v[(y*width+x)] = 0;
                }
                else {
                    _h = 0; _v = 0;
                    if (x == width - 1) {
                        h[(y*width+x)] = (-Im[(y*width+x)] + Im[(y*width)]);
                    }
                    else {
                        h[(y*width+x)] = (-Im[(y*width+x)] + Im[(y*width+x+1)]);
                    }
                    _h += pow(h[(y*width+x)], 2);

                    if (y == height - 1) {
                        v[(y*width+x)] = (-Im[(y*width+x)] + Im[x]);
                    }
                    else {
                        v[(y*width+x)] = (-Im[(y*width+x)] + Im[((y+1)*width+x)]);
                    }
                     _v += pow(v[(y*width+x)], 2);

                    if (_h + _v < lambda / beta) {
                        h[(y*width+x)] = 0;
                        v[(y*width+x)] = 0;
                    }
                }
            }
        }

        p = fftw_plan_dft_2d(height, width, Normin2In, Normin2, FFTW_FORWARD, FFTW_ESTIMATE);
        FSp = fftw_plan_dft_2d(height, width, FSIn, FS, FFTW_BACKWARD, FFTW_ESTIMATE);
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                if (x == 0) {
                    Normin2In[y*width+x][0] = h[(y*width+x)] - h[(y*width)];
                }
                else {
                    Normin2In[y*width+x][0] = h[(y*width+x-1)] - h[(y*width+x)];
                }

                if (y == 0) {
                    Normin2In[y*width+x][0] += v[(y*width+x)] - v[(x)];
                }
                else {
                    Normin2In[y*width+x][0] += v[((y-1)*width+x)] - v[(y*width+x)];
                }
            }
        }
        fftw_execute(p);
            
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                FSIn[y*width+x][0] = (Normin1[y*width+x][0] + beta*Normin2[y*width+x][0]) / Denormin[y*width+x];
                FSIn[y*width+x][1] = (Normin1[y*width+x][1] + beta*Normin2[y*width+x][1]) / Denormin[y*width+x];
            }
        }
        fftw_execute(FSp);

        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                Im[(y*width+x)] = FS[y*width+x][0] / (double)(width*height);
            }
        }

        fftw_destroy_plan(FSp);
        fftw_destroy_plan(p);

        beta = beta * kappa;
    }

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            image[y*width+x] = (unsigned char)(255. * Im[(y*width+x)]);
        }
    }

    free(h);
    free(v);
    free(Denormin);
    free(Denormin2);
    free(Ma);
    free(Im);
    fftw_free(Normin2In);
    fftw_free(FSIn);
    fftw_free(Normin1);
    fftw_free(Normin2);
    fftw_free(FS);
    
    return 0;
}

int L0smoothing_RGB(unsigned char * image, 
    int height, int width,
    unsigned char * mask, double kappa, double lambda) {

    assert(width > 0 && height > 0);
    double *Im = (double*)malloc(sizeof(double) * height * width * 3);
    for (int c = 0; c < 3; c++) {
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                Im[(y*width+x)+c*height*width] = (double)(image[3*(y*width+x) + c]) / 255.;
            }
        }
    }

    double *Ma = (double*)malloc(sizeof(double) * height * width);
    for (size_t i = 0; i < (size_t)(width * height); i++) {
        Ma[i] = (mask == NULL) ? 0. : ((double)mask[i] / 255.);
    }
    
    fftw_complex *fx_in, *fx_out;
    fftw_complex *fy_in, *fy_out;
    fftw_plan pfx, pfy;

    fx_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * width * height);
    fx_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * width * height);
    fy_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * width * height);
    fy_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * width * height);

    pfx = fftw_plan_dft_2d(height, width, fx_in, fx_out, FFTW_FORWARD, FFTW_ESTIMATE);
    pfy = fftw_plan_dft_2d(height, width, fy_in, fy_out, FFTW_FORWARD, FFTW_ESTIMATE);

    fx_in[0][0] = -1.;
    fx_in[0][1] = 0;
    fx_in[width-1][0] = 1.;
    fx_in[width-1][1] = 0;
    fy_in[0][0] = -1.;
    fy_in[0][1] = 0;
    fy_in[(height-1)*width][0] = 1.;
    fy_in[(height-1)*width][1] = 0;

    fftw_execute(pfx);
    fftw_execute(pfy);

    double *Denormin2 = (double*)malloc(sizeof(double) * width * height);
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            Denormin2[y*width + x] = pow(fx_out[y*width+x][0],2) + pow(fx_out[y*width+x][1],2)
                + pow(fy_out[y*width+x][0],2) + pow(fy_out[y*width+x][1],2);
        }
    }

    fftw_destroy_plan(pfx);
    fftw_destroy_plan(pfy);
    fftw_free(fx_in); fftw_free(fx_out);
    fftw_free(fy_in); fftw_free(fy_out);

    fftw_complex *Normin1[3], *Normin1In;
    fftw_plan p;
    Normin1In = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * width * height);
    for (int c = 0; c < 3; c++) {
        Normin1[c] = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * width * height);
        p = fftw_plan_dft_2d(height, width, Normin1In, Normin1[c], FFTW_FORWARD, FFTW_ESTIMATE);
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                Normin1In[y*width+x][0] = Im[(y*width+x)+c*height*width];
                Normin1In[y*width+x][1] = 0;
            }
        }
        fftw_execute(p);
        
        fftw_destroy_plan(p);
    }
    fftw_free(Normin1In);

    double *h = (double*)malloc(sizeof(double) * height * width * 3);
    double *v = (double*)malloc(sizeof(double) * height * width * 3);
    
    const double betamax = 1e5;
    double beta = 2 * lambda;

    double *Denormin = (double*)malloc(sizeof(double) * height * width);
    fftw_complex *FS[3], *FSIn, *Normin2[3], *Normin2In;
    fftw_plan FSp;
    FSIn = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * width * height);
    Normin2In = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * width * height);
    for (int c = 0; c < 3; c++) {
        FS[c] = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * width * height);
        Normin2[c] = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * width * height);
    }

    while (beta < betamax) {
        double _h, _v;
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                Denormin[y*width+x] = 1. + beta * Denormin2[y*width+x];

                if (Ma[(y*width+x)] > 0.5) {
                    h[(y*width+x)] = 0;
                    h[(y*width+x)+width*height] = 0;
                    h[(y*width+x)+2*width*height] = 0;
                    v[(y*width+x)] = 0;
                    v[(y*width+x)+width*height] = 0;
                    v[(y*width+x)+2*width*height] = 0;
                }
                else {
                    _h = 0; _v = 0;
                    for (int c = 0; c < 3; c++) {
                        if (x == width - 1) {
                            h[(y*width+x) + c*width*height] = (-Im[(y*width+x)+c*width*height] + Im[(y*width)+c*width*height]);
                        }
                        else {
                            h[(y*width+x) + c*width*height] = (-Im[(y*width+x)+c*width*height] + Im[(y*width+x+1)+c*width*height]);
                        }
                        _h += pow(h[(y*width+x) + c*width*height], 2);

                        if (y == height - 1) {
                            v[(y*width+x) + c*width*height] = (-Im[(y*width+x)+c*width*height] + Im[x+c*width*height]);
                        }
                        else {
                            v[(y*width+x) + c*width*height] = (-Im[(y*width+x)+c*width*height] + Im[((y+1)*width+x)+c*width*height]);
                        }
                        _v += pow(v[(y*width+x) + c*width*height], 2);
                    }

                    if (_h + _v < lambda / beta) {
                        h[(y*width+x)] = 0;
                        h[(y*width+x)+width*height] = 0;
                        h[(y*width+x)+2*width*height] = 0;
                        v[(y*width+x)] = 0;
                        v[(y*width+x)+width*height] = 0;
                        v[(y*width+x)+2*width*height] = 0;
                    }
                }
            }
        }

        for (int c = 0; c < 3; c++) {
            p = fftw_plan_dft_2d(height, width, Normin2In, Normin2[c], FFTW_FORWARD, FFTW_ESTIMATE);
            FSp = fftw_plan_dft_2d(height, width, FSIn, FS[c], FFTW_BACKWARD, FFTW_ESTIMATE);
            for (int y = 0; y < height; y++) {
                for (int x = 0; x < width; x++) {
                    if (x == 0) {
                        Normin2In[y*width+x][0] = h[(y*width+x)+c*width*height] - h[(y*width)+c*width*height];
                    }
                    else {
                        Normin2In[y*width+x][0] = h[(y*width+x-1)+c*width*height] - h[(y*width+x)+c*width*height];
                    }

                    if (y == 0) {
                        Normin2In[y*width+x][0] += v[(y*width+x)+c*width*height] - v[(x)+c*width*height];
                    }
                    else {
                        Normin2In[y*width+x][0] += v[((y-1)*width+x)+c*width*height] - v[(y*width+x)+c*width*height];
                    }
                }
            }
            fftw_execute(p);
            
            for (int y = 0; y < height; y++) {
                for (int x = 0; x < width; x++) {
                    FSIn[y*width+x][0] = (Normin1[c][y*width+x][0] + beta*Normin2[c][y*width+x][0]) / Denormin[y*width+x];
                    FSIn[y*width+x][1] = (Normin1[c][y*width+x][1] + beta*Normin2[c][y*width+x][1]) / Denormin[y*width+x];
                }
            }
            fftw_execute(FSp);

            for (int y = 0; y < height; y++) {
                for (int x = 0; x < width; x++) {
                    Im[(y*width+x)+c*width*height] = FS[c][y*width+x][0] / (double)(width*height);
                }
            }

            fftw_destroy_plan(FSp);
            fftw_destroy_plan(p);
        }

        beta = beta * kappa;
    }

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            for (int c = 0; c < 3; c++) {
                image[3*(y*width+x)+c] = (unsigned char)(255. * Im[(y*width+x) + c*width*height]);
            }
        }
    }

    free(h);
    free(v);
    free(Denormin);
    free(Denormin2);
    free(Ma);
    free(Im);
    fftw_free(Normin2In);
    fftw_free(FSIn);
    for (int c = 0; c < 3; c++) {
        fftw_free(Normin1[c]);
        fftw_free(Normin2[c]);
        fftw_free(FS[c]);
    }

    return 0;
}

int L0smoothing(unsigned char * image, 
    int height, int width, int channel,
    unsigned char * mask, double kappa, double lambda) {
    
    assert(channel == 1 || channel == 3);
    
    if (kappa == 0) {
        kappa = 2.0;
    }
    if (lambda == 0) {
        lambda = 2e-2;
    }

    if (channel == 1) {
        return L0smoothing_L(image, height, width, mask, kappa, lambda);
    }
    else if (channel == 3) {
        return L0Smoothing_RGB(image, height, width, mask, kappa, lambda);
    }
}