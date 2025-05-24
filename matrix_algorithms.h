#ifndef MATRIX_ALGORITHMS_H
#define MATRIX_ALGORITHMS_H

#include "matrix.h"

matrix *matrix_exp(const matrix *A, double eps);
int matrix_solve_gauss(matrix *A, matrix *B, matrix *X, double eps);

#endif // MATRIX_ALGORITHMS_H