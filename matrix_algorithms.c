#include "matrix_algorithms.h"
#include "matrix_ops.h"
#include "matrix_manip.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

matrix *matrix_exp(const matrix *A, double eps) {
    if (matrix_get_rows(A) != matrix_get_cols(A)) return NULL; 
    
    matrix *exp_sum = matrix_alloc_identity(matrix_get_rows(A), matrix_get_cols(A)); 
    if (!exp_sum) return NULL;
    
    matrix *current_term = matrix_copy(A); 
    if (!current_term) {
        matrix_free(exp_sum);
        return NULL;
    }
    
    int k_denominator = 1; 
    
    while (matrix_norm(current_term) >= eps) {
        matrix_add(exp_sum, current_term); 
        
        matrix *next_term_unscaled = matrix_alloc(matrix_get_rows(A), matrix_get_cols(A));
        if (!next_term_unscaled) {
            matrix_free(current_term);
            matrix_free(exp_sum);
            return NULL;
        }

        if (matrix_mul(next_term_unscaled, current_term, A) != 0) {
            matrix_free(next_term_unscaled);
            matrix_free(current_term);
            matrix_free(exp_sum);
            return NULL;
        }
        
        matrix_free(current_term); 
        current_term = next_term_unscaled; 

        k_denominator++; 

        matrix_smul(current_term, 1.0 / (double)k_denominator); 
    }
    
    matrix_free(current_term);
    return exp_sum;
}

static int find_pivot(const matrix *A, size_t column, size_t start_row, size_t *pivot_row) {
    double max_val = 0.0;
    *pivot_row = start_row;
    for (size_t i = start_row; i < matrix_get_rows(A); ++i) {
        double val = fabs(*matrix_cptr(A, i, column));
        if (val > max_val) {
            max_val = val;
            *pivot_row = i;
        }
    }
    return max_val < 1e-10 ? -1 : 0;
}

int matrix_solve_gauss(matrix *A, matrix *B, matrix *X, double eps) {
    if (matrix_get_rows(A) != matrix_get_cols(A) || matrix_get_cols(B) != 1 || matrix_get_rows(X) != matrix_get_cols(A) || matrix_get_cols(X) != 1) {
        return -1;
    }
    
    matrix *aug = matrix_alloc(matrix_get_rows(A), matrix_get_cols(A) + 1);
    if (!aug) return -1;
    for (size_t i = 0; i < matrix_get_rows(A); ++i) {
        memcpy(matrix_ptr(aug, i, 0), matrix_cptr(A, i, 0), matrix_get_cols(A) * sizeof(double));
        *matrix_ptr(aug, i, matrix_get_cols(A)) = *matrix_cptr(B, i, 0);
    }
    
    for (size_t i = 0; i < matrix_get_rows(aug); ++i) {
        size_t pivot_row;
        if (find_pivot(aug, i, i, &pivot_row) != 0) {
            matrix_free(aug);
            return -1;
        }
        if (pivot_row != i) {
            matrix_swap_rows(aug, i, pivot_row);
        }
        
        double pivot = *matrix_cptr(aug, i, i);
        if (fabs(pivot) < eps) {
            matrix_free(aug);
            return -1;
        }
        
        matrix_row_multiply(aug, i, 1.0 / pivot);
        
        for (size_t j = 0; j < matrix_get_rows(aug); ++j) {
            if (j != i) {
                double factor = *matrix_cptr(aug, j, i);
                matrix_row_add(aug, i, j, -factor);
            }
        }
    }
    
    for (size_t i = 0; i < matrix_get_rows(X); ++i) {
        *matrix_ptr(X, i, 0) = *matrix_cptr(aug, i, matrix_get_cols(A));
    }
    
    matrix_free(aug);
    return 0;
}