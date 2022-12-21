#include"guided-filter.h"

#include<stdlib.h>
#include<math.h>
#include<assert.h>
#include<limits.h>

#define util_square(a) ((a) * (a))

/*
[a.0 b.1 c.2]
[d.3 e.4 f.5]   
[g.6 h.7 i.8]
*/
void util_invMat(float M[9]) {
    static float A, B, C, D, E, F, G, H, I, det;
    A = M[4]*M[8] - M[5]*M[7];
    B = M[5]*M[6] - M[3]*M[8];
    C = M[3]*M[7] - M[4]*M[6];
    D = M[2]*M[7] - M[1]*M[8];
    E = M[0]*M[8] - M[2]*M[6];
    F = M[1]*M[6] - M[0]*M[7];
    G = M[1]*M[5] - M[2]*M[4];
    H = M[2]*M[3] - M[0]*M[5];
    I = M[0]*M[4] - M[1]*M[3];
    det = M[0] * A + M[1] * B + M[2] * C;
    M[0] = A / det; M[1] = D / det; M[2] = G / det;
    M[3] = B / det; M[4] = E / det; M[5] = H / det;
    M[6] = C / det; M[7] = F / det; M[8] = I / det;
}

int util_boxFilter(float * input,
    int height, int width, int r,
    float * output, float *tmp) {

    if (output == NULL) output = input; //overwrite input

    int tmpFlag = (tmp == NULL);
    if (tmpFlag) tmp = (float*)malloc(sizeof(float) * height * width);

    for (int x = 0; x < width; x++) {
        tmp[x] = input[x];
    }
    for (int y = 1; y < height; y++) {
        for(int x = 0; x < width; x++) {
            tmp[y*width+x] = tmp[(y-1)*width + x] + input[y*width+x];
        }
    }

    for (int y = 0; y <= r; y++) {
        for (int x = 0; x < width; x++) {
            output[y*width+x] = tmp[(y+r)*width + x];
        }
    }
    for (int y = r+1; y < height - r; y++) {
        for (int x = 0; x < width; x++) {
            output[y*width+x] = tmp[(y+r)*width + x] - tmp[(y-r-1) * width + x];
        }
    }
    for (int y = height - r; y < height; y++) {
        for (int x = 0; x < width; x++) {
            output[y*width+x] = tmp[(height-1)*width + x] - tmp[(y-r-1) * width + x];
        }
    }


    for (int y = 0; y < height; y++) {
        tmp[y*width] = output[y*width];
    }
    for (int y = 0; y < height; y++) {
        for (int x = 1; x < width; x++) {
            tmp[y*width+x] = tmp[y*width+x-1] + output[y*width+x];
        }
    }

    for (int y = 0; y < height; y++) {
        for (int x = 0; x <= r; x++) {
            output[y*width+x] = tmp[y*width+x+r];
        }
        for (int x = r+1; x < width - r; x++) {
            output[y*width+x] = tmp[y*width+x+r] - tmp[y*width+x-r-1];
        }
        for (int x = width - r; x < width; x++) {
            output[y*width+x] = tmp[y*width+width-1] - tmp[y*width+x-r-1];
        }
    }

    if (tmpFlag) free(tmp);

    return 0;
}


//overwrite p
int util_GuidedFilterMono(float * p,
    int height, int width, 
    int r, float eps,
    float * I) {

    float * tmp = (float *)malloc(sizeof(float) * width * height);

    float * N = (float *)malloc(sizeof(float) * width * height);
    for (int i = 0; i < width * height; i++) {
        N[i] = 1;
    }
    util_boxFilter(N, height, width, r, NULL, tmp);

    float * mean_I = (float *)malloc(sizeof(float) * width * height);
    float * mean_p = (float *)malloc(sizeof(float) * width * height);
    util_boxFilter(I, height, width, r, mean_I, tmp);
    util_boxFilter(p, height, width, r, mean_I, tmp);
    for (int i = 0; i < width * height; i++) {
        mean_I[i] /= N[i];
        mean_p[i] /= N[i];
    }

    //mean_Ip too
    float * cov_Ip = (float *)malloc(sizeof(float) * width * height);
    for (int i = 0; i < width * height; i++) {
        cov_Ip[i] = I[i] * p[i];
    }
    util_boxFilter(cov_Ip, height, width, r, NULL, tmp);
    for (int i = 0; i < width * height; i++) {
        cov_Ip[i] = cov_Ip[i] / N[i] - mean_I[i] * mean_p[i];
    }

    //also mean_II
    float * var_I = (float *)malloc(sizeof(float) * width * height);
    for (int i = 0; i < width * height; i++) {
        var_I[i] = util_square(I[i]);
    }
    util_boxFilter(var_I, height, width, r, NULL, tmp);
    for (int i = 0; i < width * height; i++) {
        var_I[i] = var_I[i] / N[i] - util_square(mean_I[i]);
    }

    float *mean_a = (float*)malloc(sizeof(float) * height * width);
    float *mean_b = (float*)malloc(sizeof(float) * height * width);
    for (int i = 0; i < width * height; i++) {
        mean_a[i] = cov_Ip[i] / (var_I[i] + eps);
        mean_b[i] = mean_p[i] - mean_a[i] * mean_I[i];
    }
    free(var_I);
    free(cov_Ip);
    free(mean_I); free(mean_p);
    util_boxFilter(mean_a, height, width, r, NULL, tmp);
    util_boxFilter(mean_b, height, width, r, NULL, tmp);
    for (int i = 0; i < width * height; i++) {
        mean_a[i] /= N[i];
        mean_b[i] /= N[i];
    }
    free(N);

    for (int i = 0; i < width * height; i++) {
        p[i] = mean_a[i] * I[i] + mean_b[i];
    }

    free(mean_a); free(mean_b);
    free(tmp);
    return 0;
}

int util_GuidedFilterColor(float * p,
    int height, int width, 
    int r, float eps,
    float * I) {

    float * I_r = I;
    float * I_g = I + width * height;
    float * I_b = I + 2 * width * height;

    float * tmp = (float *)malloc(sizeof(float) * width * height);

    float * N = (float *)malloc(sizeof(float) * width * height);
    for (int i = 0; i < width * height; i++) {
        N[i] = 1;
    }
    util_boxFilter(N, height, width, r, NULL, tmp);

    float * mean_p = (float *)malloc(sizeof(float) * width * height);
    util_boxFilter(p, height, width, r, mean_p, tmp);

    float * mean_I_r = (float *)malloc(sizeof(float) * width * height);
    float * mean_I_g = (float *)malloc(sizeof(float) * width * height);
    float * mean_I_b = (float *)malloc(sizeof(float) * width * height);
    util_boxFilter(I_r, height, width, r, mean_I_r, tmp);
    util_boxFilter(I_g, height, width, r, mean_I_g, tmp);
    util_boxFilter(I_b, height, width, r, mean_I_b, tmp);

    for (int i = 0; i < height * width; i++) {
        mean_p[i] /= N[i];
        mean_I_r[i] /= N[i];
        mean_I_g[i] /= N[i];
        mean_I_b[i] /= N[i];
    }

    float * cov_Ip_r = (float *)malloc(sizeof(float) * width * height);
    float * cov_Ip_g = (float *)malloc(sizeof(float) * width * height);
    float * cov_Ip_b = (float *)malloc(sizeof(float) * width * height);    
    for (int i = 0; i < width * height; i++) {
        cov_Ip_r[i] = I_r[i] * p[i];
        cov_Ip_g[i] = I_g[i] * p[i];
        cov_Ip_b[i] = I_b[i] * p[i];
    }
    util_boxFilter(cov_Ip_r, height, width, r, NULL, tmp);
    util_boxFilter(cov_Ip_g, height, width, r, NULL, tmp);
    util_boxFilter(cov_Ip_b, height, width, r, NULL, tmp);
    for (int i = 0; i < width * height; i++) {
        cov_Ip_r[i] = cov_Ip_r[i] / N[i] - mean_I_r[i] * mean_p[i];
        cov_Ip_g[i] = cov_Ip_g[i] / N[i] - mean_I_g[i] * mean_p[i];
        cov_Ip_b[i] = cov_Ip_b[i] / N[i] - mean_I_b[i] * mean_p[i];
    }

    float * var_I_rr = (float *)malloc(sizeof(float) * width * height);
    float * var_I_rg = (float *)malloc(sizeof(float) * width * height);
    float * var_I_rb = (float *)malloc(sizeof(float) * width * height);
    float * var_I_gg = (float *)malloc(sizeof(float) * width * height);
    float * var_I_gb = (float *)malloc(sizeof(float) * width * height);
    float * var_I_bb = (float *)malloc(sizeof(float) * width * height);
    for (int i = 0; i < width * height; i++) {
        var_I_rr[i] = I_r[i] * I_r[i];
        var_I_rg[i] = I_r[i] * I_g[i];
        var_I_rb[i] = I_r[i] * I_b[i];
        var_I_gg[i] = I_g[i] * I_g[i];
        var_I_gb[i] = I_g[i] * I_b[i];
        var_I_bb[i] = I_b[i] * I_b[i];
    }
    util_boxFilter(var_I_rr, height, width, r, NULL, tmp);
    util_boxFilter(var_I_rg, height, width, r, NULL, tmp);
    util_boxFilter(var_I_rb, height, width, r, NULL, tmp);
    util_boxFilter(var_I_gg, height, width, r, NULL, tmp);
    util_boxFilter(var_I_gb, height, width, r, NULL, tmp);
    util_boxFilter(var_I_bb, height, width, r, NULL, tmp);
    for (int i = 0; i < width * height; i++) {
        var_I_rr[i] = var_I_rr[i] / N[i] - mean_I_r[i] * mean_I_r[i];
        var_I_rg[i] = var_I_rg[i] / N[i] - mean_I_r[i] * mean_I_g[i];
        var_I_rb[i] = var_I_rb[i] / N[i] - mean_I_r[i] * mean_I_b[i];
        var_I_gg[i] = var_I_gg[i] / N[i] - mean_I_g[i] * mean_I_g[i];
        var_I_gb[i] = var_I_gb[i] / N[i] - mean_I_g[i] * mean_I_b[i];
        var_I_bb[i] = var_I_bb[i] / N[i] - mean_I_b[i] * mean_I_b[i];
    }

    float * a = (float*)malloc(sizeof(float) * width * height * 3);
    float Sigma[9]; //3x3 matrix
    for (int i = 0; i < width * height; i++) {
        Sigma[0] = var_I_rr[i] + eps; Sigma[1] = var_I_rg[i]      ; Sigma[2] = var_I_rb[i]      ;
        Sigma[3] = var_I_rg[i]      ; Sigma[4] = var_I_gg[i] + eps; Sigma[5] = var_I_gb[i]      ;
        Sigma[6] = var_I_rb[i]      ; Sigma[7] = var_I_gb[i]      ; Sigma[8] = var_I_bb[i] + eps;
        util_invMat(Sigma);

        a[i*3+0] = cov_Ip_r[i] * Sigma[0] + cov_Ip_g[i] * Sigma[3] + cov_Ip_b[i] * Sigma[6];
        a[i*3+1] = cov_Ip_r[i] * Sigma[1] + cov_Ip_g[i] * Sigma[4] + cov_Ip_b[i] * Sigma[7];
        a[i*3+2] = cov_Ip_r[i] * Sigma[2] + cov_Ip_g[i] * Sigma[5] + cov_Ip_b[i] * Sigma[8];
    }

    free(cov_Ip_r); free(cov_Ip_g); free(cov_Ip_b);
    free(var_I_rr); free(var_I_rg); free(var_I_rb);
    free(var_I_gg); free(var_I_gb); free(var_I_bb);

    float * b = (float*)malloc(sizeof(float) * width * height);
    for (int i = 0; i < width * height; i++) {
        b[i] = mean_p[i] - a[i*3+0] * mean_I_r[i] - a[i*3+1] * mean_I_g[i] - a[i*3+2] * mean_I_b[i];
    }
    util_boxFilter(b, height, width, r, NULL, tmp);
    
    free(mean_p);

    float *ta = (float*)malloc(sizeof(float) * width * height * 3);
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            ta[y*width+x                 ] = a[3*(y*width+x)    ] * mean_I_r[y*width+x];
            ta[y*width+x +   width*height] = a[3*(y*width+x) + 1] * mean_I_g[y*width+x];
            ta[y*width+x + 2*width*height] = a[3*(y*width+x) + 2] * mean_I_b[y*width+x];
        }
    }
    util_boxFilter(ta, height, width, r, a, tmp);
    util_boxFilter(ta + width * height, height, width, r, a + width * height, tmp);
    util_boxFilter(ta + 2 * width * height, height, width, r, a + 2 * width * height, tmp);
    
    free(ta);
    free(mean_I_r); free(mean_I_g); free(mean_I_b);

    float *a_r = a, *a_g = a + width * height, *a_b = a + 2 * width * height;

    for (int i = 0; i < width * height; i++) {
        p[i] = (a_r[i] * I_r[i] + a_g[i] * I_g[i] + a_b[i] * I_b[i] + b[i]) / N[i];
    }

    free(a); free(b);
    free(N);
    free(tmp);
    return 0;
}

int GuidedFilter(unsigned char * image,
    int height, int width, int channel,
    unsigned char * guidance,
    int guidance_channel,
    int r, float epsilon) {

    if (guidance == NULL) {
        guidance_channel = channel;
        guidance = image;
    }

    assert(height >= 10 || width >= 10);
    assert(channel == 1 || channel == 3);
    assert(guidance_channel == 1 || guidance_channel == 3);
    assert(INT_MAX / height >= width);

    float *I = NULL;
    if (guidance_channel == 1) {
        I = (float*)malloc(sizeof(float) * height * width);

        for (int i = 0; i < height * width; i++) {
            I[i] = (float)guidance[i] / 255.;
        }

    }
    else if (guidance_channel == 3) {
        I = (float*)malloc(sizeof(float) * height * width * 3);
        for (int c = 0; c < 3; c++) {
            for (int i = 0; i < height * width; i++) {
                I[i + c * width * height] = (float)guidance[3*i + c] / 255.;
            }
        }
    }

    float *p = NULL;
    if (channel == 1) {
        p = (float*)malloc(sizeof(float) * height * width);

        for (int i = 0; i < height * width; i++) {
            p[i] = (float)image[i] / 255.;
        }

        if (guidance_channel == 1) {
            util_GuidedFilterMono(p, height, width, r, epsilon, I);
        }
        else {
            util_GuidedFilterColor(p, height, width, r, epsilon, I);
        }

        float _p;
        for (int i = 0; i < width * height; i++) {
            _p = p[i];
            image[i] = (_p > 1) ? 255 : ((_p < 0) ? 0 : (unsigned char)(_p * 255));
        }
    }
    else if (channel == 3) {
        p = (float*)malloc(sizeof(float) * height * width * 3);

        for (int c = 0; c < 3; c++) {
            for (int i = 0; i < width * height; i++) {
                p[i + c*width*height] = (float)image[3*i + c] / 255.;
            }

            if (guidance_channel == 1) {
                util_GuidedFilterMono(p + c * width * height, height, width, r, epsilon, I);
            }
            else {
                util_GuidedFilterColor(p + c * width * height, height, width, r, epsilon, I);
            }
        }

        float _p;
        for (int c = 0; c < 3; c++) {
            for (int i = 0; i < width * height; i++) {
                _p = (p + c * width * height)[i];
                image[3*i + c] = (_p > 1) ? 255 : ((_p < 0) ? 0 : (unsigned char)(_p * 255));
            }
        }
    }

    free(I);
    free(p);
    return 0;
}