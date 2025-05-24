#include "matrix_manip.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

void matrix_transpose(matrix *m) {
    matrix *temp = matrix_alloc(matrix_get_cols(m), matrix_get_rows(m));
    if (!temp) return;

    for (size_t i = 0; i < matrix_get_rows(m); ++i) {
        for (size_t j = 0; j < matrix_get_cols(m); ++j) {
            *matrix_ptr(temp, j, i) = *matrix_cptr(m, i, j);
        }
    }
    matrix_replace_data(m, matrix_get_data(temp), matrix_get_rows(temp), matrix_get_cols(temp));
    matrix_free(temp);
}

#include "matrix_manip.h"
#include "matrix.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

void matrix_swap_rows(matrix *m, size_t i1, size_t i2) {
    if (i1 >= matrix_get_rows(m) || i2 >= matrix_get_rows(m)) return;
    
    const size_t cols = matrix_get_cols(m);
    double *row1 = matrix_get_data_mutable(m) + i1 * cols;
    double *row2 = matrix_get_data_mutable(m) + i2 * cols;
    
    for (size_t j = 0; j < cols; ++j) {
        double tmp = row1[j];
        row1[j] = row2[j];
        row2[j] = tmp;
    }
}

void matrix_swap_cols(matrix *m, size_t j1, size_t j2) {
    if (j1 >= matrix_get_cols(m) || j2 >= matrix_get_cols(m)) return;
    
    const size_t rows = matrix_get_rows(m);
    const size_t cols = matrix_get_cols(m);
    
    for (size_t i = 0; i < rows; ++i) {
        double *row = matrix_get_data_mutable(m) + i * cols;
        double tmp = row[j1];
        row[j1] = row[j2];
        row[j2] = tmp;
    }
}

void matrix_row_multiply(matrix *m, size_t i, double factor) {
    if (i >= matrix_get_rows(m)) return;
    
    const size_t cols = matrix_get_cols(m);
    double *row = matrix_get_data_mutable(m) + i * cols;
    
    for (size_t j = 0; j < cols; ++j) {
        row[j] *= factor;
    }
}

void matrix_row_add(matrix *m, size_t i1, size_t i2, double factor) {
    if (i1 >= matrix_get_rows(m) || i2 >= matrix_get_rows(m)) return;
    
    const size_t cols = matrix_get_cols(m);
    double *src = matrix_get_data_mutable(m) + i1 * cols;
    double *dest = matrix_get_data_mutable(m) + i2 * cols;
    
    for (size_t j = 0; j < cols; ++j) {
        dest[j] += src[j] * factor;
    }
}
double matrix_norm(const matrix *m) {
    double max_sum = 0.0;
    for (size_t i = 0; i < matrix_get_rows(m); ++i) {
        double row_sum = 0.0;
        for (size_t j = 0; j < matrix_get_cols(m); ++j) {
            row_sum += fabs(*matrix_cptr(m, i, j));
        }
        if (row_sum > max_sum) {
            max_sum = row_sum;
        }
    }
    return max_sum;
}