#ifndef MATRIX_OPS_H
#define MATRIX_OPS_H

#include "matrix.h"

int matrix_add(matrix *m1, const matrix *m2);
int matrix_sub(matrix *m1, const matrix *m2);
void matrix_smul(matrix *m, double d);
void matrix_sdiv(matrix *m, double d);

int matrix_add2(matrix *dest, const matrix *m1, const matrix *m2);
int matrix_sub2(matrix *dest, const matrix *m1, const matrix *m2);
int matrix_smul2(matrix *dest, const matrix *m, double d);
int matrix_sdiv2(matrix *dest, const matrix *m, double d);

int matrix_mul(matrix *dest, const matrix *m1, const matrix *m2);

#endif // MATRIX_OPS_H