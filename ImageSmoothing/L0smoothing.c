#include"L0smoothing.h"

#include<stdlib.h>
#include<math.h>
#include<assert.h>

#include<limits.h>

#define util_max(a, b) (((a) > (b)) ? (a) : (b))
#define util_min(a, b) (((a) < (b)) ? (a) : (b))
#define util_square(a) ((a)*(a))

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

//build bit_reverse_table among other things (if any)
int util_FFT_preload(int ** bit_reverse_table, float **root_r, float ** root_i, int n) {
    assert(n > 0);
    //assert(*bit_reverse_table == NULL); //avoid leak
    //assert(*root_r == NULL)
    //assert(*root_i == NULL)
    int size = 1;
    int s = 0, r, _s;
    while (size < n) size <<= 1, s++;

    assert(size > 0);

    *bit_reverse_table = (int*)malloc(sizeof(int) * size);
    int _i;
    for (int i = 0; i < size; i++) {
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
    *root_r = (float*)malloc(sizeof(float) * (size));
    *root_i = (float*)malloc(sizeof(float) * (size));
    float alpha, _alpha = -2 * M_PI / size;
    int j;
    for (int i = 0; i < size; i++) {
        alpha = _alpha * i;
        (*root_r)[i] = cosf(alpha);
        (*root_i)[i] = sinf(alpha);
    }

    return size;
}

int util_FFT_butterfly(float * in_r, float * in_i, int * bit_table, int n) {
    int j;
    float t;
    for (int i = 0; i < n; i++) {
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


int util_FFT_1D(float *inout_r, float *inout_i, float *root_r, float *root_i, int * bit_table, int n, int rev) {

    util_FFT_butterfly(inout_r, inout_i, bit_table, n);

    int root_idx = n;
    float _imm1_r1, _imm1_i1, _imm1_r2, _imm1_i2, real, imag, cr, ci;
    for (int block_size = 2; block_size <= n; block_size *= 2) {
        root_idx >>= 1;
        for (int block_offset = 0; block_offset < n; block_offset += block_size) {
            for (int i = 0; i < block_size / 2; i++) {
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
        for (int i = 0; i < n; i++) {
            inout_r[i] /= n;
            inout_i[i] /= n;
        }
    }

    return 0;
}

int util_FFT_preload_bluestein(float **coeff_r, float **coeff_i, int size) {
    float _alpha = (-M_PI / size), alpha;
    (*coeff_r) = (float*)malloc(sizeof(float) * size);
    (*coeff_i) = (float*)malloc(sizeof(float) * size);

    int _size = size*2;
    for (int i = 0; i < size; i++) {
        alpha = _alpha * ((i * i) % (_size));
        (*coeff_r)[i] = cosf(alpha);
        (*coeff_i)[i] = sinf(alpha);
    }

    int m = 1;
    while (m < size * 2 + 1) m<<=1;

    return m;
}

int util_FFT_1D_bluestein(float * inout_r, float * inout_i, 
    float *n_r, float * n_i,                  //coeff table
    int * bit_table, float *m_r, float *m_i,  //bit, root table
    float *temp_r, float *temp_i,             
    float *temp2_r, float *temp2_i,           
    int n, int m, int rev) {
    
    assert(n > 0 && m > 0);
    assert(m == 1 || (m & (m-1)) == 0);
    assert(m >= 2 * n + 1);

    float ar, ai, br, bi;
    
    for (int i = 0; i < n; i++) {
        ar = inout_r[i];
        ai = inout_i[i];
        br = n_r[i];
        bi = (rev) ? -n_i[i] : n_i[i];
        temp_r[i] = ar*br - ai*bi;
        temp_i[i] = ai*br + ar*bi;
    }
    for (int i = n; i < m; i++) {
        temp_r[i] = temp_i[i] = 0;
    }

    for (int i = 0; i < n; i++) {
        temp2_r[i] = n_r[i];
        temp2_i[i] = (rev) ? n_i[i] : -n_i[i];
    }
    for (int i = 0; i < m - (n * 2 - 1); i++) {
        temp2_r[i + n] = temp2_i[i + n] = 0;
    }
    for (int i = 0; i < n - 1; i++) {
        int j = i + m - n + 1;
        temp2_r[j] = n_r[n - i - 1];
        temp2_i[j] = (rev) ? n_i[n-i-1] : -n_i[n - i - 1];
    }

    util_FFT_1D(temp_r, temp_i, m_r, m_i, bit_table, m, 0);
    util_FFT_1D(temp2_r, temp2_i, m_r, m_i, bit_table, m, 0);

    for (int i = 0; i < m; i++) {
        ar = temp_r[i];
        ai = temp_i[i];
        br = temp2_r[i];
        bi = temp2_i[i];

        temp_r[i] = ar*br - ai*bi;
        temp_i[i] = ai*br + ar*bi;
    }

    util_FFT_1D(temp_r, temp_i, m_r, m_i, bit_table, m, 1);
    
    for (int i = 0; i < n; i++) {
        ar = temp_r[i];
        ai = temp_i[i];
        br = n_r[i];
        bi = (rev) ? -n_i[i] : n_i[i];
        inout_r[i] = ar*br - ai*bi;
        inout_i[i] = ai*br + ar*bi;
    }

    return 0;
}

int util_FFT_2D(float * in_r, float * in_i,
    int height, int width,
    int height_m, int width_m,
    int *height_bit_table, int *width_bit_table,
    float *height_root_r, float *height_root_i,
    float *height_coeff_r, float *height_coeff_i,
    float *height_tmp1_r, float *height_tmp1_i,
    float *height_tmp2_r, float *height_tmp2_i,
    float *width_root_r, float * width_root_i,
    float *width_coeff_r, float *width_coeff_i,
    float *width_tmp1_r, float * width_tmp1_i,
    float *width_tmp2_r, float *width_tmp2_i,
    int height_flag, int width_flag, int rev,
    float * fft_tmp_r, float *fft_tmp_i) {

    //row major to col major
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            fft_tmp_r[x*height + y] = in_r[y*width + x];
            fft_tmp_i[x*height + y] = in_i[y*width + x];
        }
    }

    //do for height
    if (height_flag) {
        for (int offset = 0; offset < width * height; offset += height) {
            util_FFT_1D(fft_tmp_r+offset, fft_tmp_i+offset, height_root_r, height_root_i, height_bit_table, height_m, rev);
        }
    }
    else {
        for (int offset = 0; offset < width * height; offset += height) {
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
        for (int offset = 0; offset < width * height; offset += width) {
            util_FFT_1D(in_r+offset, in_i+offset, width_root_r, width_root_i, width_bit_table, width_m, rev);
        }
    }
    else {
        for (int offset = 0; offset < width * height; offset += width) {
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
    unsigned char * mask, float kappa, float lambda) {

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

    float *Im = (float*)malloc(sizeof(float) * height * width * 3);
    for (int c = 0; c < 3; c++) {
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                Im[(y*width+x)+c*height*width] = (float)(image[3*(y*width+x) + c]) / 255.;
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
    float *width_root_table_r = NULL, *width_root_table_i = NULL;
    float *width_coeff_table_r = NULL, *width_coeff_table_i = NULL;
    float *width_temp1_r = NULL, *width_temp1_i = NULL;
    float *width_temp2_r = NULL, *width_temp2_i = NULL;
    int width_m;
    if (flag_width_is_power_of2) {
        width_m = width;
    }
    else {
        width_m = util_FFT_preload_bluestein(&width_coeff_table_r, &width_coeff_table_i, width);
        width_temp1_r = (float*)malloc(sizeof(float) * width_m);
        width_temp1_i = (float*)malloc(sizeof(float) * width_m);
        width_temp2_r = (float*)malloc(sizeof(float) * width_m);
        width_temp2_i = (float*)malloc(sizeof(float) * width_m);
    }
    util_FFT_preload(&width_bit_table, &width_root_table_r, &width_root_table_i, width_m);

    int flag_height_is_power_of2 = (height & (height - 1) == 0);
    int *height_bit_table = NULL;
    float *height_root_table_r = NULL, *height_root_table_i = NULL;
    float *height_coeff_table_r = NULL, *height_coeff_table_i = NULL;
    float *height_temp1_r = NULL, *height_temp1_i = NULL;
    float *height_temp2_r = NULL, *height_temp2_i = NULL;
    int height_m;
    if (flag_height_is_power_of2) {
        height_m = height;
    }
    else {
        height_m = util_FFT_preload_bluestein(&height_coeff_table_r, &height_coeff_table_i, height);
        height_temp1_r = (float*)malloc(sizeof(float) * height_m);
        height_temp1_i = (float*)malloc(sizeof(float) * height_m);
        height_temp2_r = (float*)malloc(sizeof(float) * height_m);
        height_temp2_i = (float*)malloc(sizeof(float) * height_m);
    }
    util_FFT_preload(&height_bit_table, &height_root_table_r, &height_root_table_i, height_m);

    float * fft_tmp_r = (float*)malloc(sizeof(float) * width * height); //<-- quick col-row convert
    float * fft_tmp_i = (float*)malloc(sizeof(float) * width * height);
    //

    float *fx_r = (float*)calloc(width * height, sizeof(float));
    float *fx_i = (float*)calloc(width * height, sizeof(float));
    float *fy_r = (float*)calloc(width * height, sizeof(float));
    float *fy_i = (float*)calloc(width * height, sizeof(float));
    fx_r[0] = -1;
    fx_r[width-1] = 1;
    fy_r[0] = -1;
    fy_r[(height-1) * width] = 1;

    _util_FFT2(fx_r, fx_i);
    _util_FFT2(fy_r, fy_i);

    float *Denormin2 = (float*)malloc(sizeof(float) * width * height);
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            Denormin2[y*width + x] = util_square(fx_r[y*width+x]) + util_square(fx_i[y*width+x])
                + util_square(fy_r[y*width+x]) + util_square(fy_i[y*width+x]);
        }
    }

    free(fx_r); free(fx_i);
    free(fy_r); free(fy_i);

    float *Normin1_r[3], *Normin1_i[3];
    for (int c = 0; c < 3; c++) {
        Normin1_r[c] = (float*)malloc(sizeof(float) * height * width);
        Normin1_i[c] = (float*)calloc(height * width, sizeof(float));

        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                Normin1_r[c][y * width + x] = Im[(y * width + x) + c * height * width];
            }
        }
        _util_FFT2(Normin1_r[c], Normin1_i[c]);
    }

    float *h = (float*)malloc(sizeof(float) * height * width * 3);
    float *v = (float*)malloc(sizeof(float) * height * width * 3);

    const float betamax = 1e5;
    float beta = 2 * lambda;

    //float *Denormin = (float*)malloc(sizeof(float) * height * width);
    float *FS_r = (float*)malloc(sizeof(float) * height * width);
    float *FS_i = (float*)malloc(sizeof(float) * height * width);
    float *Normin2_r = (float*)malloc(sizeof(float) * height * width);
    float *Normin2_i= (float*)malloc(sizeof(float) * height * width);

    float _Denormin;
    while (beta < betamax) {
        float _h, _v;
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
            for (int y = 0; y < height; y++) {
                for (int x = 0; x < width; x++) {
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

#undef _util_FFT2
#undef _util_IFFT2

    return 0;
}

int L0smoothing(unsigned char * image, 
    int height, int width, int channel,
    unsigned char * mask, float kappa, float lambda) {
    
    assert(width >= 10 && height >= 10);
    assert(INT_MAX / height >= width);
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