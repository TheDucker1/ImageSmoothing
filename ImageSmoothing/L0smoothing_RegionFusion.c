#include"L0smoothing_RegionFusion.h"

#include<stdlib.h>
#include<assert.h>
#include<math.h>

#include<limits.h>

#define util_max(a, b) (((a) > (b)) ? (a) : (b))
#define util_min(a, b) (((a) < (b)) ? (a) : (b))
#define util_square(a) ((a)*(a))

typedef struct _Blob Blob;
typedef struct _DoublyLinkedList DoublyLinkedList;

struct _DoublyLinkedList {
    struct _Blob * data;
    int conn;
    struct _DoublyLinkedList * next;
    struct _DoublyLinkedList * prev;
};

struct _Blob {
    int idx;
    float Y[3];
    int w;
    int * ylist;
    int * xlist;
    struct _DoublyLinkedList * nb;
};

Blob * util_create_blob(int idx, unsigned char Y[3], int y, int x) {
    Blob * b = (Blob*)malloc(sizeof(Blob));
    b->idx = idx;
    for (int i = 0; i < 3; i++) {
        b->Y[i] = (float)Y[i] / 255.f;
    }
    b->w = 1;
    b->ylist = (int*)malloc(sizeof(int));
    b->ylist[0] = y;
    b->xlist = (int*)malloc(sizeof(int));
    b->xlist[0] = x;
    b->nb = NULL;
    
    return b;
}

DoublyLinkedList* util_create_dll(Blob * data, int conn) {
    DoublyLinkedList* dll = (DoublyLinkedList*)malloc(sizeof(DoublyLinkedList));
    dll->data = data;
    dll->conn = conn;

    dll->next = dll;
    dll->prev = dll;

    return dll;
}

DoublyLinkedList* util_add_dll(DoublyLinkedList * dll, Blob * data, int conn) {
    //add when not found, otherwise add to conn
    if (dll == dll->next) {
        if (dll->data->idx == data->idx) {
            dll->conn += conn;
            return dll;
        }
        else {
            DoublyLinkedList * next = util_create_dll(data, conn);
            next->prev = dll;
            next->next = dll;
            dll->prev = next;
            dll->next = next;
            return next;
        }
    }

    while (!(dll->data->idx <= data->idx && data->idx < dll->next->data->idx)) {
        //edge case: reached the end, data->idx is either new min or max
        if (dll->data->idx > dll->next->data->idx) {
            if (data->idx > dll->data->idx || data->idx < dll->next->data->idx) {
                break;
            }
            else if (data->idx == dll->data->idx) {
                dll->conn += conn;
                return dll;
            }
            else if (data->idx == dll->next->data->idx) {
                dll->next->conn += conn;
                return dll->next;
            }
        }

        dll = dll->next;
    }

    if (dll->data->idx == data->idx) {
        dll->conn += conn;
        return dll;
    }
    else {
        DoublyLinkedList * next = util_create_dll(data, conn);
        next->prev = dll;
        next->next = dll->next;
        dll->next->prev = next;
        dll->next = next;
        return next;
    }
}

DoublyLinkedList * util_free_one_dll(DoublyLinkedList * dll) {
    if (dll == NULL) return NULL;
    if (dll == dll->next) {
        free(dll);
        return NULL;
    }

    DoublyLinkedList * next = dll->next;
    dll->prev->next = next;
    next->prev = dll->prev;
    free(dll);
    return next;
}

void util_free_all_dll(DoublyLinkedList * dll) {
    while (dll) {
        dll = util_free_one_dll(dll);
    }
}

void util_free_blob(Blob * a) {
    free(a->ylist);
    free(a->xlist);
    util_free_all_dll(a->nb);
    free(a);
}

void util_merge_blob_map(Blob * a, Blob * b) {
    int old_size = a->w;
    int new_size = a->w + b->w;

    a->ylist = realloc(a->ylist, sizeof(int) * new_size);
    a->xlist = realloc(a->xlist, sizeof(int) * new_size);

    for (int i = old_size; i < new_size; i++) {
        a->ylist[i] = b->ylist[i-old_size];
        a->xlist[i] = b->xlist[i-old_size];
    }

    a->w = new_size;
}

inline float util_cal_beta(int iter, int maxiter, float lambda) {
    return powf((float)iter / maxiter, 2.2f) * lambda;
}

int L0smoothing(unsigned char * image, 
    int height, int width, int channel,
    unsigned char * mask, float lambda) {

    assert(width >= 10 && height >= 10);
    assert(channel = 3);
    assert(INT_MAX / height >= width);
    assert(lambda > 0);

    //initial
    int M = height * width;
    int P = M;
    
    Blob ** blob_map = (Blob**)malloc(sizeof(Blob*) * M);
    Blob *blob_cur, *blob_tmp;
    int idx = 0;
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            blob_cur = blob_map[idx] = util_create_blob(idx, image+idx*3, y, x);
            if (x > 0) {
                blob_tmp = blob_map[idx - 1];
                if (blob_cur->nb == NULL) {
                    blob_cur->nb = util_create_dll(blob_tmp, 1);
                }
                else {
                    blob_cur->nb = util_add_dll(blob_cur->nb, blob_tmp, 1);
                }

                if (blob_tmp->nb == NULL) {
                    blob_tmp->nb = util_create_dll(blob_cur, 1);
                }
                else {
                    blob_tmp->nb = util_add_dll(blob_tmp->nb, blob_cur, 1);
                }
            }

            if (y > 0) {
                blob_tmp = blob_map[idx - width];
                if (blob_cur->nb == NULL) {
                    blob_cur->nb = util_create_dll(blob_tmp, 1);
                }
                else {
                    blob_cur->nb = util_add_dll(blob_cur->nb, blob_tmp, 1);
                }

                if (blob_tmp->nb == NULL) {
                    blob_tmp->nb = util_create_dll(blob_cur, 1);
                }
                else {
                    blob_tmp->nb = util_add_dll(blob_tmp->nb, blob_cur, 1);
                }
            }

            idx++;
        }
    }

    int iter = 0, maxiter = 50, test = 0;
    float beta = util_cal_beta(iter, maxiter, lambda);
    while (beta <= lambda && P > 1) { //beta will > lambda when iter > maxiter

        for (int i = 0; i < M; i++) {
            if (blob_map[i] == NULL) {
                continue;
            }

            blob_cur = blob_map[i];
            DoublyLinkedList * nb_anchor = blob_cur->nb;

            do {

                blob_tmp = blob_cur->nb->data;
                int wi = blob_cur->w, wj = blob_tmp->w;
                int c = blob_cur->nb->conn;
                double dif = 0;
                for (int j = 0; j < 3; j++) {
                    dif += util_square(blob_cur->Y[j] - blob_tmp->Y[j]);
                }
                if (dif * wi * wj <= c * beta * (wi + wj)) {

                    for (int j = 0; j < 3; j++) {
                        blob_cur->Y[j] = (wi * blob_cur->Y[j] + wj * blob_tmp->Y[j]) / (wi + wj);
                    }

                    util_merge_blob_map(blob_cur, blob_tmp);

                    DoublyLinkedList * cur_nb = blob_cur->nb;
                    DoublyLinkedList * tmp_nb = blob_tmp->nb;

                    //fomd tmp in cur_nb (should not move)
                    //cur_nb = util_add_dll(cur_nb, blob_tmp, 0);
                    cur_nb = util_free_one_dll(cur_nb);
                    
                    //find cur in tmp_nb
                    tmp_nb = util_add_dll(tmp_nb, blob_cur, 0);
                    tmp_nb = util_free_one_dll(tmp_nb);
                    DoublyLinkedList *anchor = tmp_nb;
                    if (tmp_nb) {
                        int cjk = 0;
                        do {
                            //add new nb to cur
                            if (cur_nb) {
                                cur_nb = util_add_dll(cur_nb, tmp_nb->data, tmp_nb->conn);
                            }
                            else {
                                cur_nb = util_create_dll(tmp_nb->data, tmp_nb->conn);
                            }

                            //remove tmp from nb
                            tmp_nb->data->nb = util_add_dll(tmp_nb->data->nb, blob_tmp, 0);
                            cjk = tmp_nb->conn;
                            tmp_nb->data->nb = util_free_one_dll(tmp_nb->data->nb);

                            //add cur to new nb
                            if (tmp_nb->data->nb) {
                                tmp_nb->data->nb = util_add_dll(tmp_nb->data->nb, blob_cur, cjk);
                            }
                            else {
                                tmp_nb->data->nb = util_create_dll(blob_cur, cjk);
                            }
                            tmp_nb = tmp_nb->next;
                        } while (tmp_nb != anchor);
                    }
                    
                    blob_cur->nb = cur_nb;
                    blob_tmp->nb = tmp_nb;

                    //remove blob_tmp
                    blob_map[blob_tmp->idx] = NULL;
                    util_free_blob(blob_tmp);

                    P--;
                    break;
                }

                blob_cur->nb = blob_cur->nb->next;
            } while (blob_cur->nb != nb_anchor);
        }

        iter++;
        beta = util_cal_beta(iter, maxiter, lambda);
    }

    int check = 0;
    for (int i = 0; i < M; i++) {
        if (blob_map[i]) {
            blob_cur = blob_map[i];
            check += blob_cur->w;
            unsigned char c[3];
            float _c;
            for (int k = 0; k < 3; k++) {
                _c = blob_cur->Y[k];
                c[k] = (_c < 0) ? 0 : ((_c > 1) ? 255 : (unsigned char)(_c * 255));
            }
            for (int j = 0; j < blob_cur->w; j++) {
                int idx = blob_cur->ylist[j] * width + blob_cur->xlist[j];
                for (int k = 0; k < 3; k++) {
                    image[3*idx + k] = c[k];
                }
            }

            util_free_blob(blob_map[i]);
        }
    }

    free(blob_map);

    return 0;
}