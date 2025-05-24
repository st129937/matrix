#include "matrix_ops.h"
#include "matrix.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

static int check_same_size(const matrix *m1, const matrix *m2) {
    return (matrix_get_rows(m1) == matrix_get_rows(m2)) && 
           (matrix_get_cols(m1) == matrix_get_cols(m2));
}

static int check_mult_size(const matrix *m1, const matrix *m2) {
    return matrix_get_cols(m1) == matrix_get_rows(m2);
}

int matrix_add(matrix *m1, const matrix *m2) {
    if (!check_same_size(m1, m2)) return -1;
    
    const size_t elements = matrix_get_rows(m1) * matrix_get_cols(m1);
    double *data1 = matrix_get_data_mutable(m1);
    const double *data2 = matrix_get_data(m2);
    
    for (size_t i = 0; i < elements; ++i) {
        data1[i] += data2[i];
    }
    return 0;
}

int matrix_sub(matrix *m1, const matrix *m2) {
    if (!check_same_size(m1, m2)) return -1;
    
    const size_t elements = matrix_get_rows(m1) * matrix_get_cols(m1);
    double *data1 = matrix_get_data_mutable(m1);
    const double *data2 = matrix_get_data(m2);
    
    for (size_t i = 0; i < elements; ++i) {
        data1[i] -= data2[i];
    }
    return 0;
}

void matrix_smul(matrix *m, double d) {
    const size_t elements = matrix_get_rows(m) * matrix_get_cols(m);
    double *data = matrix_get_data_mutable(m);
    
    for (size_t i = 0; i < elements; ++i) {
        data[i] *= d;
    }
}

void matrix_sdiv(matrix *m, double d) {
    matrix_smul(m, 1.0 / d);
}

int matrix_add2(matrix *dest, const matrix *m1, const matrix *m2) {
    if (!check_same_size(m1, m2) || 
        !check_same_size(m1, dest)) return -1;
    
    const size_t elements = matrix_get_rows(m1) * matrix_get_cols(m1);
    const double *data1 = matrix_get_data(m1);
    const double *data2 = matrix_get_data(m2);
    double *dest_data = matrix_get_data_mutable(dest);
    
    for (size_t i = 0; i < elements; ++i) {
        dest_data[i] = data1[i] + data2[i];
    }
    return 0;
}

int matrix_sub2(matrix *dest, const matrix *m1, const matrix *m2) {
    if (!check_same_size(m1, m2) || 
        !check_same_size(m1, dest)) return -1;
    
    const size_t elements = matrix_get_rows(m1) * matrix_get_cols(m1);
    const double *data1 = matrix_get_data(m1);
    const double *data2 = matrix_get_data(m2);
    double *dest_data = matrix_get_data_mutable(dest);
    
    for (size_t i = 0; i < elements; ++i) {
        dest_data[i] = data1[i] - data2[i];
    }
    return 0;
}

int matrix_smul2(matrix *dest, const matrix *m, double d) {
    if (!check_same_size(m, dest)) return -1;
    
    const size_t elements = matrix_get_rows(m) * matrix_get_cols(m);
    const double *src_data = matrix_get_data(m);
    double *dest_data = matrix_get_data_mutable(dest);
    
    for (size_t i = 0; i < elements; ++i) {
        dest_data[i] = src_data[i] * d;
    }
    return 0;
}

int matrix_sdiv2(matrix *dest, const matrix *m, double d) {
    return matrix_smul2(dest, m, 1.0 / d);
}

int matrix_mul(matrix *dest, const matrix *m1, const matrix *m2) {
    if (!check_mult_size(m1, m2)) return -1;
    if (matrix_get_rows(dest) != matrix_get_rows(m1) || 
        matrix_get_cols(dest) != matrix_get_cols(m2)) return -1;
    
    const size_t m = matrix_get_rows(m1);
    const size_t n = matrix_get_cols(m2);
    const size_t p = matrix_get_cols(m1);
    
    matrix *temp = matrix_alloc(m, n);
    if (!temp) return -1;
    
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            double sum = 0.0;
            for (size_t k = 0; k < p; ++k) {
                sum += *matrix_cptr(m1, i, k) * *matrix_cptr(m2, k, j);
            }
            *matrix_ptr(temp, i, j) = sum;
        }
    }
    
    int result = matrix_assign(dest, temp);
    matrix_free(temp);
    return result;
}