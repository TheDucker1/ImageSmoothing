#include"eap.h"

#include<stdlib.h>
#include<math.h>
#include<assert.h>
#include<stdint.h>
#include<time.h>

#define eapSquare(a) ((a)*(a))

/*-----https://stackoverflow.com/a/26237777------*/
static unsigned int g_seed;

// Used to seed the generator.           
inline void fast_srand(int seed) {
    g_seed = seed;
}

// Compute a pseudorandom integer.
// Output value in range [0, 32767]
inline double fast_rand(void) {
    g_seed = (214013*g_seed+2531011);
    return (double)((g_seed>>16)&0x7FFF + 1) / 32769;
}
/*-----------------------------------------------*/

int cmpFunc(const void * a, const void * b) {
    return *(double*)a > *(double*)b;
}

inline void quickPartitionSwap(double * arr, int i, int j, double t) {
    t = arr[i];
    arr[i] = arr[j];
    arr[j] = t;
}

void quickPartition(double * arr, int lo, int hi, int idx) {
    if (lo >= hi) {
        return;
    }
    double t=0;
    int left = lo-1, right = hi;
    //pivot = arr[hi]
    for (int i = lo; i < right; i++) {
        if (arr[i] < arr[hi]) {
            left++;
            quickPartitionSwap(arr, i, left, t);
        }

        else if (arr[i] > arr[hi]) {
            right--;
            quickPartitionSwap(arr, i, right, t);
            i--;
        }
    }
    quickPartitionSwap(arr, hi, right, t);

    if (left < idx && idx <= right) {
        return;
    }
    else if (idx <= left) {
        return quickPartition(arr, lo, left, idx);
    }
    else {
        return quickPartition(arr, right+1, hi, idx);
    }
}

int eap(unsigned char * image, 
    int height, int width,
    double u,
    int (*func)(unsigned char *, int, int, unsigned char*)) {

    assert(height > 0 && width > 0);
    assert(0. <= u && u <= 1.);

    const int iter = 5;

    const double eps_weight = 1.0;

    fast_srand(time(NULL));

    unsigned char *output = (unsigned char*)malloc(sizeof(unsigned char) * height * width * 3);
    unsigned char *Ma = (unsigned char*)malloc(sizeof(unsigned char) * height * width);
    double *value = (double*)malloc(sizeof(double) * height * width);
    double *weight = (double*)malloc(sizeof(double) * height * width);
    double *thre = (double*)malloc(sizeof(double) * height * width);
    double r, alpha, _knapsack, _thre, r2;

    int ny, nx;

    for (size_t i = 0; i < (size_t)(width) * height; i++) {
        Ma[i] = (fast_rand() < 0.5) ? 0 : 255;
    }

    for (int k = 1; k <= iter; k++) {

        for (size_t i = 0; i < (size_t)(width) * height * 3; i++) {
            output[i] = image[i];
        }

        func(output, height, width, Ma);
        
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                value[y*width+x] = 0;
                weight[y*width+x] = 0;
                for (int c = 0; c < 3; c++) {
                    r =  0;
                    
                    for (int dy = -1; dy <= 1; dy++) {
                        for (int dx = -1; dx <= 1; dx++) {
                            ny = y + dy;
                            if (ny == -1) {
                                ny = 0;
                            }
                            else if (ny == height) {
                                ny = height-1;
                            }

                            nx = x + dx;
                            if (nx == -1) {
                                nx = 0;
                            }
                            else if (nx == width) {
                                nx = width - 1;
                            }

                            r += output[3*(ny*width+nx) + c];
                        }
                    }
                    
                    r /= 9;
                    r2 = (double)output[3*(y*width+x)+c] / 255. - image[3*(y*width+x)+c]; //don't know why but the matlab implementation imply it takes the euclid distance of [0, 1] out_img and [0, 255] in_img??
                    r -= output[3*(y*width+x)+c];

                    value[y*width+x] += eapSquare(r2);
                    weight[y*width+x] += eapSquare(r);
                }

                value[y*width+x] /= (weight[y*width+x] + eps_weight);
                thre[y*width+x] = value[y*width+x];
            }
        }

        int idx = (int)((double)(height*width) * (1. - (double)k * u / iter));
        
        quickPartition(thre, 0, height*width-1, idx);
        //qsort(thre, height * width, sizeof(double), cmpFunc);
        _thre = thre[idx];

        alpha = ((double)(iter) - k) / (iter-1);
        for (size_t i = 0; i < (size_t)height*width; i++) {
            _knapsack = value[i] / _thre;

            Ma[i] = (fast_rand() * alpha + 1. - alpha < _knapsack) ? 255 : 0;
        }
    }

    for (size_t i = 0; i < height * width * 3; i++) {
        image[i] = output[i];
    }

    free(output);
    free(Ma);
    free(value);
    free(weight);
    free(thre);

    return 0;
}