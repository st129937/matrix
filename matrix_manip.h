#ifndef MATRIX_MANIP_H
#define MATRIX_MANIP_H

#include "matrix.h"

void matrix_transpose(matrix *m);
void matrix_swap_rows(matrix *m, size_t i1, size_t i2);
void matrix_swap_cols(matrix *m, size_t j1, size_t j2);
void matrix_row_multiply(matrix *m, size_t i, double factor);
void matrix_row_add(matrix *m, size_t i1, size_t i2, double factor);
double matrix_norm(const matrix *m);

#endif // MATRIX_MANIP_H