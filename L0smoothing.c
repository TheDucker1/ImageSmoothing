#include"L0smoothing.h"

#include<stdlib.h>
#include<math.h>
#include<assert.h>

//#include"fftw3.h" //a more practical version

#define util_max(a, b) (((a) > (b)) ? (a) : (b))
#define util_min(a, b) (((a) < (b)) ? (a) : (b))
#define util_square(a) ((a)*(a))

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define _util_FFT2(a, b) util_FFT_2D(a, b, height, width, height_m, width_m, \
        height_bit_table, width_bit_table, \
        height_root_table_r, height_root_table_i, \
        height_coeff_table_r, height_coeff_table_i, \
        height_temp1_r, height_temp1_i, \
        height_temp2_r, height_temp2_i, \
        width_root_table_r, width_root_table_i, \
        width_coeff_table_r, width_coeff_table_i, \
        width_temp1_r, width_temp1_i, \
        width_temp2_r, width_temp2_i, \
        flag_height_is_power_of2, flag_width_is_power_of2, \
        0, \
        fft_tmp_r, fft_tmp_i)

#define _util_IFFT2(a, b) util_FFT_2D(a, b, height, width, height_m, width_m, \
        height_bit_table, width_bit_table, \
        height_root_table_r, height_root_table_i, \
        height_coeff_table_r, height_coeff_table_i, \
        height_temp1_r, height_temp1_i, \
        height_temp2_r, height_temp2_i, \
        width_root_table_r, width_root_table_i, \
        width_coeff_table_r, width_coeff_table_i, \
        width_temp1_r, width_temp1_i, \
        width_temp2_r, width_temp2_i, \
        flag_height_is_power_of2, flag_width_is_power_of2, \
        1, \
        fft_tmp_r, fft_tmp_i)

//TO DO: simple functions to replace fftw

//build bit_reverse_table among other things (if any)
long util_FFT_preload(int ** bit_reverse_table, double **root_r, double ** root_i, long n) {
    assert(n > 0);
    //assert(*bit_reverse_table == NULL); //avoid leak
    //assert(*root_r == NULL)
    //assert(*root_i == NULL)
    long size = 1;
    int s = 0, r, _s;
    while (size < n) size <<= 1, s++;

    assert(size > 0);

    *bit_reverse_table = (int*)malloc(sizeof(int) * size);
    long _i;
    for (long i = 0; i < size; i++) {
        r = 0;
        _i = i;
        _s = s;
        for (; _i; _i >>= 1) {
            r <<= 1;
            r |= _i & 1;
            _s--;
        }
        r <<= _s;
        (*bit_reverse_table)[i] = r;
    }

    //pre calculate root
    *root_r = (double*)malloc(sizeof(double) * (size));
    *root_i = (double*)malloc(sizeof(double) * (size));
    double alpha, _alpha = -2 * M_PI / size;
    long j;
    for (long i = 0; i < size; i++) {
        alpha = _alpha * i;
        (*root_r)[i] = cos(alpha);
        (*root_i)[i] = sin(alpha);
    }

    return size;
}

int util_FFT_butterfly(double * in_r, double * in_i, int * bit_table, long n) {
    long j;
    double t;
    for (long i = 0; i < n; i++) {
        j = bit_table[i];
        if (i < j) {
            t = in_r[i];
            in_r[i] = in_r[j];
            in_r[j] = t;

            t = in_i[i];
            in_i[i] = in_i[j];
            in_i[j] = t;
        }
    }

    return 0;
}


int util_FFT_1D(double *inout_r, double *inout_i, double *root_r, double *root_i, int * bit_table, long n, int rev) {

    util_FFT_butterfly(inout_r, inout_i, bit_table, n);

    long root_idx = n;
    double _imm1_r1, _imm1_i1, _imm1_r2, _imm1_i2, real, imag, cr, ci;
    for (long block_size = 2; block_size <= n; block_size *= 2) {
        root_idx >>= 1;
        for (long block_offset = 0; block_offset < n; block_offset += block_size) {
            for (long i = 0; i < block_size / 2; i++) {
                _imm1_r1 = inout_r[block_offset + i];
                _imm1_i1 = inout_i[block_offset + i];
                _imm1_r2 = inout_r[block_offset + i + block_size/2];
                _imm1_i2 = inout_i[block_offset + i + block_size/2];
                cr = root_r[i * root_idx];
                ci = (rev) ? -root_i[i * root_idx] : root_i[i * root_idx];
                real = _imm1_r2 * cr - _imm1_i2 * ci;
                imag = _imm1_r2 * ci + _imm1_i2 * cr;

                inout_r[block_offset + i] = _imm1_r1 + real;
                inout_i[block_offset + i] = _imm1_i1 + imag;
                inout_r[block_offset + i + block_size/2] = _imm1_r1 - real;
                inout_i[block_offset + i + block_size/2] = _imm1_i1 - imag;
            }
        }
    }

    if (rev) {
        for (long i = 0; i < n; i++) {
            inout_r[i] /= n;
            inout_i[i] /= n;
        }
    }

    return 0;
}

long util_FFT_preload_bluestein(double **coeff_r, double **coeff_i, long size) {
    double _alpha = (-M_PI / size), alpha;
    (*coeff_r) = (double*)malloc(sizeof(double) * size);
    (*coeff_i) = (double*)malloc(sizeof(double) * size);

    long _size = size*2;
    for (long i = 0; i < size; i++) {
        alpha = _alpha * ((i * i) % (_size));
        (*coeff_r)[i] = cos(alpha);
        (*coeff_i)[i] = sin(alpha);
    }

    long m = 1;
    while (m < size * 2 + 1) m<<=1;

    return m;
}

int util_FFT_1D_bluestein(double * inout_r, double * inout_i, 
    double *n_r, double * n_i,                  //coeff table
    int * bit_table, double *m_r, double *m_i,  //bit, root table
    double *temp_r, double *temp_i,             
    double *temp2_r, double *temp2_i,           
    long n, long m, int rev) {
    
    assert(n > 0 && m > 0);
    assert(m == 1 || (m & (m-1)) == 0);
    assert(m >= 2 * n + 1);

    double ar, ai, br, bi;
    
    for (long i = 0; i < n; i++) {
        ar = inout_r[i];
        ai = inout_i[i];
        br = n_r[i];
        bi = (rev) ? -n_i[i] : n_i[i];
        temp_r[i] = ar*br - ai*bi;
        temp_i[i] = ai*br + ar*bi;
    }
    for (long i = n; i < m; i++) {
        temp_r[i] = temp_i[i] = 0;
    }

    for (long i = 0; i < n; i++) {
        temp2_r[i] = n_r[i];
        temp2_i[i] = (rev) ? n_i[i] : -n_i[i];
    }
    for (long i = 0; i < m - (n * 2 - 1); i++) {
        temp2_r[i + n] = temp2_i[i + n] = 0;
    }
    for (long i = 0; i < n - 1; i++) {
        long j = i + m - n + 1;
        temp2_r[j] = n_r[n - i - 1];
        temp2_i[j] = (rev) ? n_i[n-i-1] : -n_i[n - i - 1];
    }

    util_FFT_1D(temp_r, temp_i, m_r, m_i, bit_table, m, 0);
    util_FFT_1D(temp2_r, temp2_i, m_r, m_i, bit_table, m, 0);

    for (long i = 0; i < m; i++) {
        ar = temp_r[i];
        ai = temp_i[i];
        br = temp2_r[i];
        bi = temp2_i[i];

        temp_r[i] = ar*br - ai*bi;
        temp_i[i] = ai*br + ar*bi;
    }

    util_FFT_1D(temp_r, temp_i, m_r, m_i, bit_table, m, 1);
    
    for (long i = 0; i < n; i++) {
        ar = temp_r[i];
        ai = temp_i[i];
        br = n_r[i];
        bi = (rev) ? -n_i[i] : n_i[i];
        inout_r[i] = ar*br - ai*bi;
        inout_i[i] = ai*br + ar*bi;
    }

    return 0;
}

int util_FFT_2D(double * in_r, double * in_i,
    long height, long width,
    long height_m, long width_m,
    int *height_bit_table, int *width_bit_table,
    double *height_root_r, double *height_root_i,
    double *height_coeff_r, double *height_coeff_i,
    double *height_tmp1_r, double *height_tmp1_i,
    double *height_tmp2_r, double *height_tmp2_i,
    double *width_root_r, double * width_root_i,
    double *width_coeff_r, double *width_coeff_i,
    double *width_tmp1_r, double * width_tmp1_i,
    double *width_tmp2_r, double *width_tmp2_i,
    int height_flag, int width_flag, int rev,
    double * fft_tmp_r, double *fft_tmp_i) {

    //row major to col major
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            fft_tmp_r[x*height + y] = in_r[y*width + x];
            fft_tmp_i[x*height + y] = in_i[y*width + x];
        }
    }

    //do for height
    if (height_flag) {
        for (long offset = 0; offset < width * height; offset += height) {
            util_FFT_1D(fft_tmp_r+offset, fft_tmp_i+offset, height_root_r, height_root_i, height_bit_table, height_m, rev);
        }
    }
    else {
        for (long offset = 0; offset < width * height; offset += height) {
            util_FFT_1D_bluestein(fft_tmp_r+offset, fft_tmp_i+offset, 
                height_coeff_r, height_coeff_i, height_bit_table, 
                height_root_r, height_root_i,
                height_tmp1_r, height_tmp1_i,
                height_tmp2_r, height_tmp2_i,
                height, height_m, rev);
        }
    }

    //col major back to row major
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            in_r[y*width + x] = fft_tmp_r[x*height + y];
            in_i[y*width + x] = fft_tmp_i[x*height + y];
        }
    }

    //do for width
    if (width_flag) {
        for (long offset = 0; offset < width * height; offset += width) {
            util_FFT_1D(in_r+offset, in_i+offset, width_root_r, width_root_i, width_bit_table, width_m, rev);
        }
    }
    else {
        for (long offset = 0; offset < width * height; offset += width) {
            util_FFT_1D_bluestein(in_r+offset, in_i+offset, 
                width_coeff_r, width_coeff_i, width_bit_table, 
                width_root_r, width_root_i,
                width_tmp1_r, width_tmp1_i,
                width_tmp2_r, width_tmp2_i,
                width, width_m, rev);
        }
    }


    return 0;
}

int L0smoothing_RGB(unsigned char * image, 
    int height, int width,
    unsigned char * mask, double kappa, double lambda) {

    double *Im = (double*)malloc(sizeof(double) * height * width * 3);
    for (int c = 0; c < 3; c++) {
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                Im[(y*width+x)+c*height*width] = (double)(image[3*(y*width+x) + c]) / 255.;
            }
        }
    }

    unsigned char *Ma = (unsigned char*)malloc(sizeof(unsigned char) * height * width);
    for (size_t i = 0; i < (size_t)(width * height); i++) {
        Ma[i] = (mask == NULL) ? 0 : mask[i];
    }

    //FFT routine
    int flag_width_is_power_of2 = (width & (width - 1) == 0);
    int *width_bit_table = NULL;
    double *width_root_table_r = NULL, *width_root_table_i = NULL;
    double *width_coeff_table_r = NULL, *width_coeff_table_i = NULL;
    double *width_temp1_r = NULL, *width_temp1_i = NULL;
    double *width_temp2_r = NULL, *width_temp2_i = NULL;
    long width_m;
    if (flag_width_is_power_of2) {
        width_m = width;
    }
    else {
        width_m = util_FFT_preload_bluestein(&width_coeff_table_r, &width_coeff_table_i, width);
        width_temp1_r = (double*)malloc(sizeof(double) * width_m);
        width_temp1_i = (double*)malloc(sizeof(double) * width_m);
        width_temp2_r = (double*)malloc(sizeof(double) * width_m);
        width_temp2_i = (double*)malloc(sizeof(double) * width_m);
    }
    util_FFT_preload(&width_bit_table, &width_root_table_r, &width_root_table_i, width_m);

    int flag_height_is_power_of2 = (height & (height - 1) == 0);
    int *height_bit_table = NULL;
    double *height_root_table_r = NULL, *height_root_table_i = NULL;
    double *height_coeff_table_r = NULL, *height_coeff_table_i = NULL;
    double *height_temp1_r = NULL, *height_temp1_i = NULL;
    double *height_temp2_r = NULL, *height_temp2_i = NULL;
    long height_m;
    if (flag_height_is_power_of2) {
        height_m = height;
    }
    else {
        height_m = util_FFT_preload_bluestein(&height_coeff_table_r, &height_coeff_table_i, height);
        height_temp1_r = (double*)malloc(sizeof(double) * height_m);
        height_temp1_i = (double*)malloc(sizeof(double) * height_m);
        height_temp2_r = (double*)malloc(sizeof(double) * height_m);
        height_temp2_i = (double*)malloc(sizeof(double) * height_m);
    }
    util_FFT_preload(&height_bit_table, &height_root_table_r, &height_root_table_i, height_m);

    double * fft_tmp_r = (double*)malloc(sizeof(double) * width * height); //<-- quick col-row convert
    double * fft_tmp_i = (double*)malloc(sizeof(double) * width * height);
    //

    double *fx_r = (double*)calloc(width * height, sizeof(double));
    double *fx_i = (double*)calloc(width * height, sizeof(double));
    double *fy_r = (double*)calloc(width * height, sizeof(double));
    double *fy_i = (double*)calloc(width * height, sizeof(double));
    fx_r[0] = -1;
    fx_r[width-1] = 1;
    fy_r[0] = -1;
    fy_r[(height-1) * width] = 1;

    _util_FFT2(fx_r, fx_i);
    _util_FFT2(fy_r, fy_i);

    double *Denormin2 = (double*)malloc(sizeof(double) * width * height);
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            Denormin2[y*width + x] = util_square(fx_r[y*width+x]) + util_square(fx_i[y*width+x])
                + util_square(fy_r[y*width+x]) + util_square(fy_i[y*width+x]);
        }
    }

    free(fx_r); free(fx_i);
    free(fy_r); free(fy_i);

    double *Normin1_r[3], *Normin1_i[3];
    for (int c = 0; c < 3; c++) {
        Normin1_r[c] = (double*)malloc(sizeof(double) * height * width);
        Normin1_i[c] = (double*)calloc(height * width, sizeof(double));

        for (long y = 0; y < height; y++) {
            for (long x = 0; x < width; x++) {
                Normin1_r[c][y * width + x] = Im[(y * width + x) + c * height * width];
            }
        }
        _util_FFT2(Normin1_r[c], Normin1_i[c]);
    }

    double *h = (double*)malloc(sizeof(double) * height * width * 3);
    double *v = (double*)malloc(sizeof(double) * height * width * 3);

    const double betamax = 1e5;
    double beta = 2 * lambda;

    //double *Denormin = (double*)malloc(sizeof(double) * height * width);
    double *FS_r = (double*)malloc(sizeof(double) * height * width);
    double *FS_i = (double*)malloc(sizeof(double) * height * width);
    double *Normin2_r = (double*)malloc(sizeof(double) * height * width);
    double *Normin2_i= (double*)malloc(sizeof(double) * height * width);

    double _Denormin;
    while (beta < betamax) {
        double _h, _v;
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                //Denormin[y*width+x] = 1. + beta * Denormin2[y*width+x];

                if (Ma[(y*width+x)] > 127) {
                    h[(y*width+x)] = h[(y*width+x)+width*height] = h[(y*width+x)+2*width*height] = 0;
                    v[(y*width+x)] = v[(y*width+x)+width*height] = v[(y*width+x)+2*width*height] = 0;
                }
                else {
                    _h = _v = 0;
                    for (int c = 0; c < 3; c++) {
                        if (x == width - 1) {
                            h[(y*width+x) + c*width*height] = (-Im[(y*width+x)+c*width*height] + Im[(y*width)+c*width*height]);
                        }
                        else {
                            h[(y*width+x) + c*width*height] = (-Im[(y*width+x)+c*width*height] + Im[(y*width+x+1)+c*width*height]);
                        }
                        _h += util_square(h[(y*width+x) + c*width*height]);

                        if (y == height - 1) {
                            v[(y*width+x) + c*width*height] = (-Im[(y*width+x)+c*width*height] + Im[x+c*width*height]);
                        }
                        else {
                            v[(y*width+x) + c*width*height] = (-Im[(y*width+x)+c*width*height] + Im[((y+1)*width+x)+c*width*height]);
                        }
                        _v += util_square(v[(y*width+x) + c*width*height]);
                    }

                    if (_h + _v < lambda / beta) {
                        h[(y*width+x)] = h[(y*width+x)+width*height] = h[(y*width+x)+2*width*height] = 0;
                        v[(y*width+x)] = v[(y*width+x)+width*height] = v[(y*width+x)+2*width*height] = 0;
                    }
                }
            }
        }

        for (int c = 0; c < 3; c++) {
            for (long y = 0; y < height; y++) {
                for (long x = 0; x < width; x++) {
                    Normin2_i[y*width + x] = 0;
                    if (x == 0) {
                        Normin2_r[y*width + x] = h[(y*width+width-1)+c*width*height] - h[(y*width)+c*width*height];
                    }
                    else {
                        Normin2_r[y*width + x] = h[(y*width+x-1)+c*width*height] - h[(y*width+x)+c*width*height];
                    }

                    if (y == 0) {
                        Normin2_r[y*width + x] += v[((height-1)*width+x)+c*width*height] - v[(x)+c*width*height];
                    }
                    else {
                        Normin2_r[y*width + x] += v[((y-1)*width+x)+c*width*height] - v[(y*width+x)+c*width*height];
                    }
                }
            }

            _util_FFT2(Normin2_r, Normin2_i);

            
            for (int y = 0; y < height; y++) {
                for (int x = 0; x < width; x++) {
                    _Denormin = 1. + beta * Denormin2[y*width+x];
                    FS_r[y*width+x] = (Normin1_r[c][y*width+x] + beta*Normin2_r[y*width+x]) / _Denormin;
                    FS_i[y*width+x] = (Normin1_i[c][y*width+x] + beta*Normin2_i[y*width+x]) / _Denormin;
                }
            }

            _util_IFFT2(FS_r, FS_i);

            for (int y = 0; y < height; y++) {
                for (int x = 0; x < width; x++) {
                    Im[(y*width+x)+c*width*height] = FS_r[y*width+x] / (width * height);
                }
            }
        }

        beta *= kappa;
    }

    free(FS_r); free(FS_i);
    free(Normin2_r); free(Normin2_i);

    //free(Denormin);

    for (int c = 0; c < 3; c++) {
        free(Normin1_r[c]); free(Normin1_i[c]);
    }

    free(fft_tmp_r);
    free(fft_tmp_i);

    free(width_bit_table);
    free(width_root_table_r);
    free(width_root_table_i);
    free(width_coeff_table_r);
    free(width_coeff_table_i);
    free(width_temp1_r);
    free(width_temp1_i);
    free(width_temp2_r);
    free(width_temp2_i);
    
    free(height_bit_table);
    free(height_root_table_r);
    free(height_root_table_i);
    free(height_coeff_table_r);
    free(height_coeff_table_i);
    free(height_temp1_r);
    free(height_temp1_i);
    free(height_temp2_r);
    free(height_temp2_i);

    free(Ma);

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            for (int c = 0; c < 3; c++) {
                image[3*(y*width+x)+c] = (Im[(y*width+x) + c*width*height] > 1) ? 255 : ((Im[(y*width+x) + c*width*height] < 0) ? 0 : (unsigned char)(255 * Im[(y*width+x) + c*width*height]));
            }
        }
    }
    free(Im);

    return 0;
}

/* //require FFTW
int _L0smoothing_RGB(unsigned char * image, 
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

    unsigned char *Ma = (unsigned char*)malloc(sizeof(unsigned char) * height * width);
    for (size_t i = 0; i < (size_t)(width * height); i++) {
        Ma[i] = (mask == NULL) ? 0 : mask[i];
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

    for (size_t i = 0; i < width * height; i++) {
        fx_in[i][0] = fx_in[i][1] = fy_in[i][0] = fy_in[i][1] = 0;
    }

    fx_in[0][0] = -1.;
    fx_in[width-1][0] = 1.;
    fy_in[0][0] = -1.;
    fy_in[(height-1)*width][0] = 1.;

    fftw_execute(pfx);
    fftw_execute(pfy);

    double *Denormin2 = (double*)malloc(sizeof(double) * width * height);
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            Denormin2[y*width + x] = util_square(fx_out[y*width+x][0]) + util_square(fx_out[y*width+x][1])
                + util_square(fy_out[y*width+x][0]) + util_square(fy_out[y*width+x][1]);
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

                if (Ma[(y*width+x)] > 127) {
                    h[(y*width+x)] = 0;
                    h[(y*width+x)+width*height] = 0;
                    h[(y*width+x)+2*width*height] = 0;
                    v[(y*width+x)] = 0;
                    v[(y*width+x)+width*height] = 0;
                    v[(y*width+x)+2*width*height] = 0;
                }
                else {
                    _h = _v = 0;
                    for (int c = 0; c < 3; c++) {
                        if (x == width - 1) {
                            h[(y*width+x) + c*width*height] = (-Im[(y*width+x)+c*width*height] + Im[(y*width)+c*width*height]);
                        }
                        else {
                            h[(y*width+x) + c*width*height] = (-Im[(y*width+x)+c*width*height] + Im[(y*width+x+1)+c*width*height]);
                        }
                        _h += util_square(h[(y*width+x) + c*width*height]);

                        if (y == height - 1) {
                            v[(y*width+x) + c*width*height] = (-Im[(y*width+x)+c*width*height] + Im[x+c*width*height]);
                        }
                        else {
                            v[(y*width+x) + c*width*height] = (-Im[(y*width+x)+c*width*height] + Im[((y+1)*width+x)+c*width*height]);
                        }
                        _v += util_square(v[(y*width+x) + c*width*height]);
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
                    Normin2In[y*width+x][1] = 0;
                    if (x == 0) {
                        Normin2In[y*width+x][0] = h[(y*width+width-1)+c*width*height] - h[(y*width)+c*width*height];
                    }
                    else {
                        Normin2In[y*width+x][0] = h[(y*width+x-1)+c*width*height] - h[(y*width+x)+c*width*height];
                    }

                    if (y == 0) {
                        Normin2In[y*width+x][0] += v[((height-1)*width+x)+c*width*height] - v[(x)+c*width*height];
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
                    Im[(y*width+x)+c*width*height] = FS[c][y*width+x][0] / (width*height);
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
                image[3*(y*width+x)+c] = (Im[(y*width+x) + c*width*height] > 1) ? 255 : ((Im[(y*width+x) + c*width*height] < 0) ? 0 : (unsigned char)(255 * Im[(y*width+x) + c*width*height]));
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
*/

int L0smoothing(unsigned char * image, 
    int height, int width, int channel,
    unsigned char * mask, double kappa, double lambda) {
    
    assert(width >= 10 && height >= 10);
    assert(channel == 3);

    if (kappa == 0) {
        kappa = 2.0;
    }
    if (lambda == 0) {
        lambda = 2e-2;
    }
    
    if (channel == 3) {
        return L0smoothing_RGB(image, height, width, mask, kappa, lambda);
    }

    return 0;
}