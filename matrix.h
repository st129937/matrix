#ifndef MATRIX_H
#define MATRIX_H

#include <stddef.h>

typedef struct matrix matrix;

const double* matrix_get_data(const matrix *m);
size_t matrix_get_rows(const matrix *m);
size_t matrix_get_cols(const matrix *m);
double* matrix_get_data_mutable(matrix *m);

int matrix_replace_data(matrix *m, const double *new_data, size_t rows, size_t cols);

matrix *matrix_alloc(size_t rows, size_t cols);
matrix *matrix_copy(const matrix *m);
void matrix_free(matrix *m);

double *matrix_ptr(matrix *m, size_t i, size_t j);
const double *matrix_cptr(const matrix *m, size_t i, size_t j);

void matrix_set_zero(matrix *m);
void matrix_set_identity(matrix *m);

matrix *matrix_alloc_zero(size_t rows, size_t cols);
matrix *matrix_alloc_identity(size_t rows, size_t cols);

int matrix_assign(matrix *dest, const matrix *src);

void matrix_print(const matrix *m);
int matrix_read(matrix *m);

#endif // MATRIX_H