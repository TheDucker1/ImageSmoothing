#include"eap.h"

#include<stdlib.h>
#include<math.h>
#include<assert.h>
#include<stdint.h>
#include<time.h>

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
    return (double)((g_seed>>16)&0x7FFF) / 32767;
}
/*-----------------------------------------------*/

void quickPartitionSwap(double * arr, size_t i, size_t j) {
    double t;
    t = arr[i];
    arr[i] = arr[j];
    arr[j] = t;
}

void quickPartition(double * arr, size_t lo, size_t hi, size_t idx) {
    if (lo >= hi) {
        return;
    }
    size_t left = lo-1, right = hi;
    //pivot = arr[hi]
    for (size_t i = lo; i < right; i++) {
        if (arr[i] < arr[hi]) {
            left++;
            quickPartitionSwap(arr, i, left);
        }

        else if (arr[i] > arr[hi]) {
            right--;
            quickPartitionSwap(arr, i, right);
            i--;
        }
    }
    quickPartitionSwap(arr, hi, right);

    if (left < idx && idx <= right) {
        return;
    }
    else if (idx <= left) {
        quickPartition(arr, lo, left, idx);
    }
    else {
        quickPartition(arr, right+1, hi, idx);
    }
}

int eap(unsigned char * image, 
    int height, int width,
    double u, int iter,
    int (*func)(unsigned char *, int, int, unsigned char*)) {

    assert(height > 0 && width > 0);
    assert(0. <= u && u <= 1.);

    if (iter == 0) {
        iter = 5;
    }

    const double eps_weight = 1.0;

    fast_srand(time(NULL));

    unsigned char *output = (unsigned char*)malloc(sizeof(unsigned char) * height * width * 3);
    unsigned char *Ma = (unsigned char*)malloc(sizeof(unsigned char) * height * width);
    double *value = (double*)malloc(sizeof(double) * height * width);
    double *weight = (double*)malloc(sizeof(double) * height * width);
    double *knapsack = (double*)malloc(sizeof(double) * height * width);
    double *thre = (double*)malloc(sizeof(double) * height * width);
    double r, alpha, _knapsack, _thre;

    int ny, nx;

    for (size_t i = 0; i < (size_t)(width) * height; i++) {
        r = fast_rand();
        Ma[i] = (r < 0.5) ? 0 : 255;
    }
     
    for (int k = 0; k < iter; k++) {
        
        for (size_t i = 0; i < (size_t)(width) * height * 3; i++) {
            output[i] = image[i];
        }

        func(output, height, width, Ma);
        
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                value[y*width+x] = 0;
                weight[y*width+x] = 0;
                for (int c = 0; c < 3; c++) {
                    r = 0;
                    
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

                            r += (double)output[3*(ny*width+nx) + c];
                        }
                    }
                    
                    r /= 9;
                    value[y*width+x] += pow(((double)(output[3*(y*width+x) + c]) - image[3*(y*width+x) + c]), 2);
                    weight[y*width+x] += pow(((double)(output[3*(y*width+x) + c]) - r), 2);
                }

                knapsack[y*width+x] = value[y*width+x] / (weight[y*width+x] + eps_weight);
                thre[y*width+x] = knapsack[y*width+x];
            }
        }

        alpha = ((double)(iter) - (k+1)) / (iter-1);

        size_t idx = (size_t)((double)(height*width) * (1. - ((double)(k+1) * u) / (iter)));
        quickPartition(thre, 0, height*width-1, idx);
        _thre = thre[idx];
        
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                _knapsack = knapsack[y*width+x] / _thre;
                r = fast_rand();

                if (r * alpha + 1. - alpha < _knapsack) {
                    Ma[y*width+x] = 255;
                }
                else {
                    Ma[y*width+x] = 0;
                }
            }
        }

    }

    for (size_t i = 0; i < height * width * 3; i++) {
        image[i] = output[i];
    }

    free(output);
    free(Ma);
    free(value);
    free(weight);
    free(knapsack);
    free(thre);

    return 0;
}