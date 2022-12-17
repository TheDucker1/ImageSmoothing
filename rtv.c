#include"rtv.h"

#include<stdlib.h>
#include<math.h>
#include<assert.h>

#include<limits.h>

#define util_max(a, b) (((a) > (b)) ? (a) : (b))
#define util_min(a, b) (((a) < (b)) ? (a) : (b))
#define util_square(a) ((a)*(a))
#define util_inv(a) (1. / (a))  //imply double

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int rtv_RGB_A_mul(double * A_D, double * A_dy, double * A_dx, int height,
    int k, double * in, double * out) {

    // out = A * in

    // dx; dy; D; -dy; -dx
    int f[4];
    for (int i = 0; i < k; i++) {
        f[0] = i >= height;
        f[1] = i >= 1;
        f[2] = i + 1 < k;
        f[3] = i + height < k;

        out[i] = A_D[i] * in[i];
        for (int j = 0; j < 4; j++) {
            if (f[j]) {
                switch (j) {
                    case 0:
                    out[i] += A_dx[i-height] * in[i-height];
                    break;

                    case 1:
                    out[i] += A_dy[i-1] * in[i-1];
                    break;

                    case 2:
                    out[i] += A_dy[i] * in[i+1];
                    break;

                    case 3:
                    out[i] += A_dx[i] * in[i+height];
                    break;
                }
            }
        }
    }

    return 0;
}

int rtv_RGB_precond(double * L_D, double * L_dy, double * L_dx, int height,
    int k, double * b, double * x) {

    // M = L * L'
    // M * x = b

    double _x;
    for (int i = 0; i < k; i++) {
        _x = b[i];
        if (i >= 1) _x -= x[i-1] * L_dy[i-1];
        if (i >= height) _x -= x[i-height] * L_dx[i-height];

        x[i] = _x / L_D[i];
    }

    for (int i = k; i >= 0; i--) {
        _x = x[i];

        if (i < k - 1) _x -= L_dy[i] * x[i+1];
        if (i < k - height) _x -= L_dx[i] * x[i+height];

        x[i] = _x / L_D[i];
    }

    return 0;
}

int rtv_RGB_solveLinearEquation(double * IN,
    int height, int width,
    double * wx, double * wy,
    double lambda,
    double * OUT) {
    
    int k = height * width;

    double * D = (double*)malloc(sizeof(double) * k);
    double * dx = (double*)malloc(sizeof(double) * k);
    double * dy = (double*)malloc(sizeof(double) * k);

    double * OUT_p;
    double * IN_p;

    double _e, _s, _w, _n; // 4dir
    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            _e = -lambda * wx[x*height + y];
            _s = -lambda * wy[x*height + y];
            _w = (x != 0) ? (dx[(x-1) * height + y]) : 0;
            _n = (x == 0 && y == 0) ? 0 : dy[x * height + y - 1];

            dx[x*height + y] = _e;
            dy[x*height + y] = _s;
            D[x*height + y] = 1 - (_e + _w + _s + _n);
        }
    }

    /*
        D on main diag
    A = dx on -r, +r diag
        dy on -1, +1 diag
    */
    // dx[0:r]
    // dy[0:-1]
    // D
    // dy[0:-1]
    // dx[0:r]

    //CALCULATE start
    //IC start
    double * L_D = (double*)malloc(sizeof(double) * k);
    double * L_dy = (double*)malloc(sizeof(double) * (k - 1));
    double * L_dx = (double*)malloc(sizeof(double) * (k - height));

    double _L_D;
    for (size_t i = 0; i < k; i++) {
        _L_D = D[i];
        if (i >= 1) _L_D -= util_square(L_dy[i-1]);
        if (i >= height) _L_D -= util_square(L_dx[i-height]);

        _L_D = sqrt(_L_D);
        L_D[i] = _L_D;

        //assuming height > 2
        if (i < k - 1) {
            L_dy[i] = dy[i] / _L_D;
            //L_dy[i] = (dy[i] - dx[i] * dy[i-1]) / _L_D (height == 2, should not occur)
        }

        if (i < k - height) {
            L_dx[i] = dx[i] / _L_D;
        }
    }
    //IC end

    //PCG start
    //https://www.math.uci.edu/~chenint/CAMtips/CG.html
    //function u = pcg(A,b,u,B,tol)
    // u: OUT_p
    // b: IN_p
    double normb, normr, rho, rho_old, pAp, tol = 1e-1, alpha, beta;
    double * r = (double *)malloc(sizeof(double) * k);
    double * p = (double *)malloc(sizeof(double) * k);
    double * Ap = (double *)malloc(sizeof(double) * k);
    double * Br = (double *)malloc(sizeof(double) * k);
    int step = 0, maxstep = 100;

    for (int c = 0; c < 3; c++) {
        IN_p = IN + c * width * height;
        OUT_p = OUT + c * width * height;

        step = 1;

        //initial guess
        normb = normr = 0;
        for (int i = 0; i < k; i++) {
            normb += util_square(IN_p[i]);
            OUT_p[i] = 0;
            r[i] = IN_p[i];
        }
        normr = normb = sqrt(normb);
        normb *= tol;

        //start computing
        rho = 1e10;
        while (normr > normb && step <= maxstep) {
            rtv_RGB_precond(L_D, L_dy, L_dx, height, k, r, Br);
            
            rho = 0;
            for (int i = 0; i < k; i++) {
                rho += r[i] * (Br[i]);
            }

            if (step == 1) {
                for (int i = 0; i < k; i++) {
                    p[i] = Br[i];
                }
            }
            else {
                beta = rho / rho_old;
                for (int i = 0; i < k; i++) {
                    p[i] = Br[i] + beta * p[i];
                }
            }

            rtv_RGB_A_mul(D, dy, dx, height, k, p, Ap);
            pAp = 0;
            for (int i = 0; i < k; i++) {
                pAp += p[i] * (Ap[i]);
            }

            alpha = rho / pAp;

            normr = 0;
            for (int i = 0; i < k; i++) {
                OUT_p[i] += alpha * p[i];
                r[i] -= alpha * Ap[i];
                normr += util_square(r[i]);
            }
            normr = sqrt(normr);

            rho_old = rho;
            step++;
        }
    }

    free(r);
    free(p);
    free(Ap);
    free(Br);
    //PCG end
    //CALCULATE end

    free(D);
    free(dx);
    free(dy);
    free(L_D);
    free(L_dx);
    free(L_dy);


    return 0;
}

int rtv_RGB_computeTextureWeights(double * fin,
    int height, int width,
    double sigma, double sharpness,
    double *retx, double * rety) {
    
    double vareps = 0.001, _fx, _fy, _wto = 0, _fbin = 0;

    double *wto = (double*)malloc(sizeof(double) * width * height);
    double *fbin = (double*)malloc(sizeof(double) * width * height * 3);
    double *fbin2 = (double*)malloc(sizeof(double) * width * height * 3);

    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            _wto = 0;
            for (int c = 0; c < 3; c++) {
                _fx = (x == width-1) ? 0 : (fin[((x+1)*height+y) + c*width*height] - fin[(x*height+y) + c*width*height]);
                _fy = (y == height-1) ? 0 : (fin[(x*height+y+1) + c*width*height] - fin[(x*height+y) + c*width*height]);
                _wto += sqrt(util_square(_fx) + util_square(_fy));
            }
            _wto = util_max(_wto / 3, sharpness);

            wto[x*height+y] = util_inv(_wto);
        }
    }


    //START OF lpfilter
    int ksize = (int)round(5*sigma) | 1;
    int r = ksize / 2;
    double * kernel = (double*)malloc(sizeof(double) * ksize), s = 0;
    for (int i = 0; i <= r; i++) {
        kernel[r-i] = kernel[r+i] = util_inv(exp((double)(i*i) / (2 * (sigma*sigma)))) * util_inv(sigma * sqrt(2 * M_PI));
        if (i == 0) s += kernel[r-i];
        else s += 2 * kernel[r-i];
    }
    
    for (int i = 0; i < ksize; i++) {
        kernel[i] /= s;
    }
    

    for (int c = 0; c < 3; c++) {
        for (int x = 0; x < width; x++) { //first pass, y
            for (int y = 0; y < height; y++) {
                _fbin = 0;
                int ny;
                for (int dy = -r; dy <= r; dy++) {
                    ny = y + dy;
                    if (ny < 0) ny = -ny - 1;
                    if (ny >= height) {
                        ny = ((ny / height) % 2) ? height - (ny % height) - 1 : ny % height;
                    }
                    //if (ny < 0) ny = 0;
                    //else if (ny >= height) ny = height-1;

                    _fbin += fin[x*height + ny + c*width*height] * kernel[dy+r];
                }
                fbin2[x*height + y + c*width*height] = _fbin;
            }
        }

        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) { //second pass, x
                _fbin = 0;
                int nx;
                for (int dx = -r; dx <= r; dx++) {
                    nx = x + dx;
                    if (nx < 0) nx = -nx - 1;
                    if (nx >= width) {
                        nx = ((nx / width) % 2) ? width - (nx % width) - 1 : nx % width;
                    }
                    //if (nx < 0) nx = 0;
                    //else if (nx >= width) nx = width - 1;

                    _fbin += fbin2[nx*height+y + c*width*height] * kernel[dx+r];
                }
                fbin[x*height+y + c *width*height] = _fbin;
            }
        }
    }

    
    free(fbin2);
    free(kernel);

    //END OF lpfilter

    double * fbin_pointer2[3];
    fbin_pointer2[0] = fbin;
    fbin_pointer2[1] = fbin + width * height;
    fbin_pointer2[2] = fbin + 2 * width * height;

    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            //use _fx and _fy
            _fx = _fy = 0;
            
            if (x != width - 1) {
                for (int c = 0; c < 3; c++) {
                    _fx += fabs(fbin_pointer2[c][x*height+y] - fbin_pointer2[c][(x+1)*height+y]);
                }
            }

            if (y != height - 1) {
                for (int c = 0; c < 3; c++) {
                    _fy += fabs(fbin_pointer2[c][x*height+y] - fbin_pointer2[c][x*height+y+1]);
                }
            }

            _fx = util_inv(util_max(_fx / 3, vareps));
            _fy = util_inv(util_max(_fy / 3, vareps));

            retx[x*height+y] = (x == width - 1) ? 0 : _fx * wto[x*height+y];
            rety[x*height+y] = (y == height - 1) ? 0 : _fy * wto[x*height+y];

        }
    }

    free(fbin);
    free(wto);

    return 0;
}

int rtv_RGB(unsigned char * image, 
    int height, int width,
    unsigned char * mask,
    double lambda, double sigma, double sharpness, int maxIter) {

    double * I = (double*)malloc(sizeof(double) * width * height * 3);
    double * wx = (double*)malloc(sizeof(double) * width * height);
    double * wy = (double*)malloc(sizeof(double) * width * height);
    double * Ix = (double*)malloc(sizeof(double) * width * height * 3);

    if (mask == NULL) {
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                for (int c = 0; c < 3; c++) {
                    I[3*(y*width+x) + c] = (double)image[3*(y*width+x) + c];
                }
            }
        }
    }
    else {
        double * M = (double*)malloc(sizeof(double) * width * height);
        double * Mw = (double*)malloc(sizeof(double) * width * height);

        double _Mw, _M;

        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                M[y*width+x] = (double)mask[y*width + x] / 255.;
                Mw[y*width+x] = 1.0 - M[y*width+x] + 1e-5;
                for (int c = 0; c < 3; c++) {
                    I[3*(y*width+x) + c] = (double)image[3*(y*width+x) + c];
                }
            }
        }

        int piter = 32;
        for (int k = 1; k <= piter; k++) {
            for (int y = 0; y < height; y++) {
                for (int x = 0; x < width; x++) {
                    _Mw = Mw[y*width+x];
                    _M = M[y*width+x];
                    for (int c = 0; c < 3; c++) {
                        I[3*(y*width+x) + c] += I[3*(y*width+x) + c];
                        if (x == width - 1) I[3*(y*width+x) + c] += I[3*(y*width+x) + c];
                        else I[3*(y*width+x) + c] += I[3*(y*width+x+1) + c];
                        if (y == height - 1) I[3*(y*width+x) + c] += I[3*(y*width+x) + c];
                        else I[3*(y*width+x) + c] += I[3*((y+1)*width+x) + c];

                        I[3*(y*width+x) + c] *= _Mw;
                    }

                    _Mw += Mw[y*width+x];
                    if (x == width - 1) _Mw += Mw[y*width+x];
                    else _Mw += Mw[y*width+x+1];
                    if (y == height - 1) _Mw += Mw[y*width+x];
                    else _Mw += Mw[(y+1)*width+x];

                    for (int c = 0; c < 3; c++) {
                        I[3*(y*width+x) + c] /= _Mw;
                        I[3*(y*width+x) + c] = I[3*(y*width+x) + c] * _M + (double)image[3*(y*width+x) + c] * (1. - _M);
                    }

                    Mw[y*width+x] = (_Mw / 4 < 0.25) ? 1e-5 : (1 + 1e-5);
                }
            }
        }
        free(M);
        free(Mw);
    }

    for (size_t i = 0; i < (size_t)height * width * 3; i++) {
        if (I[i] > 255) {
            I[i] = 1; 
        }
        else if (I[i] < 0) {
            I[i] = 0;
        }
        else {
            I[i] /= 255;
        }
    }

    double sigma_iter = sigma;
    lambda /= 2;
    double dec = 2;

    //transform to col major, channel by channel
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            for (int c = 0; c < 3; c++) {
                Ix[(x*height+y) + c * width * height] = I[3*(y*width+x) + c];
            }
        }
    }

    for (size_t i = 0; i < height * width * 3; i++) {
        I[i] = Ix[i];
    }
    //I have the form (channel1 | channel2 | channel3)

    for (int iter = 1; iter <= maxIter; iter++) {
        rtv_RGB_computeTextureWeights(Ix, height, width, sigma_iter, sharpness, wx, wy);
        rtv_RGB_solveLinearEquation(I, height, width, wx, wy, lambda, Ix);

        sigma_iter /= dec;
        if (sigma_iter < 0.5) sigma_iter = 0.5;
    }

    double _I;
    for (int c = 0; c < 3; c++) {
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                _I = Ix[(x*height+y) + c*width * height];
                image[3*(y*width+x) + c] = (_I > 1) ? 255 : ((_I < 0) ? 0 : (unsigned char)(_I * 255));
            }
        }
    }

    free(I);
    free(Ix);
    free(wx);
    free(wy);

    return 0;
}

int rtv(unsigned char * image, 
    int height, int width, int channel,
    unsigned char * mask,
    double lambda, double sigma, double sharpness, int maxIter) {

    assert(channel == 3);
    assert(width >= 10 && height >= 10);
    assert(INT_MAX / height > width);

    if (lambda == 0) lambda = 0.01;
    if (sigma == 0) sigma = 3;
    if (sharpness == 0) sharpness = 0.02;
    if (maxIter == 0) maxIter = 4;

    return rtv_RGB(image, height, width, mask, lambda, sigma, sharpness, maxIter);

}