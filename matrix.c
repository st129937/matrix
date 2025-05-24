#include "matrix.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

struct matrix {
    double *data;
    size_t rows;
    size_t cols;
};

const double* matrix_get_data(const matrix *m) {
    return m->data;
}

double* matrix_get_data_mutable(matrix *m) {
    return m->data;
}

size_t matrix_get_rows(const matrix *m) {
    return m->rows;
}

size_t matrix_get_cols(const matrix *m) {
    return m->cols;
}

int matrix_replace_data(matrix *m, const double *new_data, size_t rows, size_t cols) {
    if (rows * cols != m->rows * m->cols) return -1;
    
    memcpy(m->data, new_data, rows * cols * sizeof(double));
    m->rows = rows;
    m->cols = cols;
    return 0;
}

matrix *matrix_alloc(size_t rows, size_t cols) {
    if (rows == 0 || cols == 0) return NULL;
    matrix *m = malloc(sizeof(matrix));
    if (!m) return NULL;
    m->rows = rows;
    m->cols = cols;
    m->data = malloc(rows * cols * sizeof(double));
    if (!m->data) {
        free(m);
        return NULL;
    }
    return m;
}

matrix *matrix_copy(const matrix *m) {
    if (!m) return NULL;
    matrix *copy = matrix_alloc(m->rows, m->cols);
    if (!copy) return NULL;
    memcpy(copy->data, m->data, m->rows * m->cols * sizeof(double));
    return copy;
}

void matrix_free(matrix *m) {
    if (m) {
        free(m->data);
        free(m);
    }
}

double *matrix_ptr(matrix *m, size_t i, size_t j) {
    return &m->data[i * m->cols + j];
}

const double *matrix_cptr(const matrix *m, size_t i, size_t j) {
    return &m->data[i * m->cols + j];
}

void matrix_set_zero(matrix *m) {
    memset(m->data, 0, m->rows * m->cols * sizeof(double));
}

void matrix_set_identity(matrix *m) {
    matrix_set_zero(m);
    size_t min_dim = m->rows < m->cols ? m->rows : m->cols;
    for (size_t i = 0; i < min_dim; ++i) {
        *matrix_ptr(m, i, i) = 1.0;
    }
}

matrix *matrix_alloc_zero(size_t rows, size_t cols) {
    matrix *m = matrix_alloc(rows, cols);
    if (m) matrix_set_zero(m);
    return m;
}

matrix *matrix_alloc_identity(size_t rows, size_t cols) {
    matrix *m = matrix_alloc(rows, cols);
    if (m) matrix_set_identity(m);
    return m;
}

int matrix_assign(matrix *dest, const matrix *src) {
    if (!dest || !src) return -1;
    if (dest->rows != src->rows || dest->cols != src->cols) return -1;
    memcpy(dest->data, src->data, dest->rows * dest->cols * sizeof(double));
    return 0;
}

void matrix_print(const matrix *m) {
    for (size_t i = 0; i < m->rows; ++i) {
        for (size_t j = 0; j < m->cols; ++j) {
            printf("%8.4f ", *matrix_cptr(m, i, j));
        }
        printf("\n");
    }
}

int matrix_read(matrix *m) {
    for (size_t i = 0; i < m->rows; ++i) {
        for (size_t j = 0; j < m->cols; ++j) {
            if (scanf("%lf", matrix_ptr(m, i, j)) != 1) {
                return -1;
            }
        }
    }
    return 0;
}