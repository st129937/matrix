#include "matrix.h"
#include "matrix_ops.h"
#include "matrix_manip.h"
#include "matrix_algorithms.h"
#include <stdio.h>
#include <locale.h>

void print_menu() {
    printf("\n=== Меню операций с матрицами ===\n");
    printf("1. Показать исходные матрицы\n");
    printf("2. Сложить матрицы\n");
    printf("3. Вычесть матрицы\n");
    printf("4. Умножить на скаляр\n");
    printf("5. Транспонировать матрицу\n");
    printf("6. Поменять строки местами\n");
    printf("7. Поменять столбцы местами\n");
    printf("8. Умножить матрицы\n");
    printf("9. Скопировать матрицу\n");
    printf("10. Матричная экспонента\n");
    printf("11. Решить СЛАУ Ax=B (тест Гаусса)\n");
    printf("0. Выход\n");
    printf("Выберите операцию: ");
}

void print_matrix_choice() {
    printf("\nВыберите матрицу:\n");
    printf("1. Матрица A\n");
    printf("2. Матрица B\n");
    printf("Ваш выбор: ");
}

int main() {
    setlocale(LC_CTYPE, "Russian");
    
    matrix *A = matrix_alloc(2, 2);
    *matrix_ptr(A, 0, 0) = 2; *matrix_ptr(A, 0, 1) = 0;
    *matrix_ptr(A, 1, 0) = 0; *matrix_ptr(A, 1, 1) = -1;

    matrix *B = matrix_alloc(2, 2);
    *matrix_ptr(B, 0, 0) = 5; *matrix_ptr(B, 0, 1) = 6;
    *matrix_ptr(B, 1, 0) = 7; *matrix_ptr(B, 1, 1) = 8;

    matrix *result = NULL;
    int choice, m_choice;
    double scalar;

    do {
        print_menu();
        scanf("%d", &choice);
        
        switch(choice) {
            case 1:
                printf("\nМатрица A:\n");
                matrix_print(A);
                printf("\nМатрица B:\n");
                matrix_print(B);
                break;
                
            case 2: 
                result = matrix_copy(A);
                if(matrix_add(result, B)) {
                    printf("Ошибка: разные размеры матриц!\n");
                    matrix_free(result);
                } else {
                    printf("\nA + B =\n");
                    matrix_print(result);
                    matrix_free(result);
                }
                break;
                
            case 3: 
                result = matrix_copy(A);
                if(matrix_sub(result, B)) {
                    printf("Ошибка: разные размеры матриц!\n");
                    matrix_free(result);
                } else {
                    printf("\nA - B =\n");
                    matrix_print(result);
                    matrix_free(result);
                }
                break;
                
            case 4: 
                printf("Введите скаляр: ");
                scanf("%lf", &scalar);
                result = matrix_copy(A);
                matrix_smul(result, scalar);
                printf("\nA * %.2f =\n", scalar);
                matrix_print(result);
                matrix_free(result);
                break;
                
            case 5: 
                result = matrix_copy(A);
                matrix_transpose(result);
                printf("\nA^T =\n");
                matrix_print(result);
                matrix_free(result);
                break;
                
            case 6: 
                result = matrix_copy(A);
                matrix_swap_rows(result, 0, 1);
                printf("\nПосле перестановки строк:\n");
                matrix_print(result);
                matrix_free(result);
                break;
                
            case 7: 
                result = matrix_copy(A);
                matrix_swap_cols(result, 0, 1);
                printf("\nПосле перестановки столбцов:\n");
                matrix_print(result);
                matrix_free(result);
                break;
                
            case 8: 
                result = matrix_alloc(2, 2);
                if(matrix_mul(result, A, B)) {
                    printf("Ошибка умножения!\n");
                } else {
                    printf("\nA * B =\n");
                    matrix_print(result);
                }
                matrix_free(result);
                break;
                
            case 9:
                print_matrix_choice();
                scanf("%d", &m_choice);
                
                if(m_choice == 1) {
                    result = matrix_copy(A);
                    printf("\nКопия матрицы A:\n");
                } else if(m_choice == 2) {
                    result = matrix_copy(B);
                    printf("\nКопия матрицы B:\n");
                } else {
                    result = NULL;
                    printf("Неверный выбор!\n");
                }
                
                matrix_print(result);
                matrix_free(result);
                break;
                
            case 10: 
                result = matrix_exp(A, 1e-6);
                if(result) {
                    printf("\nexp(A) =\n");
                    matrix_print(result);
                    matrix_free(result);
                } else {
                    printf("Ошибка вычисления экспоненты!\n");
                }
                break;
            
            case 11:
                {
                    matrix *gauss_A = matrix_alloc(3, 3);
                    matrix *gauss_B = matrix_alloc(3, 1);
                    matrix *gauss_X = matrix_alloc(3, 1);

                    if (!gauss_A || !gauss_B || !gauss_X) {
                        if(gauss_A) matrix_free(gauss_A);
                        if(gauss_B) matrix_free(gauss_B);
                        if(gauss_X) matrix_free(gauss_X);
                        break;
                    }

                    *matrix_ptr(gauss_A, 0, 0) = 0; *matrix_ptr(gauss_A, 0, 1) = 0; *matrix_ptr(gauss_A, 0, 2) = 1;
                    *matrix_ptr(gauss_A, 1, 0) = 0; *matrix_ptr(gauss_A, 1, 1) = 1; *matrix_ptr(gauss_A, 1, 2) = 0;
                    *matrix_ptr(gauss_A, 2, 0) = 1; *matrix_ptr(gauss_A, 2, 1) = 0; *matrix_ptr(gauss_A, 2, 2) = 0;

                    *matrix_ptr(gauss_B, 0, 0) = 1;
                    *matrix_ptr(gauss_B, 1, 0) = 2;
                    *matrix_ptr(gauss_B, 2, 0) = 3;

                    printf("Матрица A для теста:\n");
                    matrix_print(gauss_A);
                    printf("Вектор B для теста:\n");
                    matrix_print(gauss_B);

                    matrix_solve_gauss(gauss_A, gauss_B, gauss_X, 1e-9);
                    printf("Решение X:\n");
                    matrix_print(gauss_X);

                    matrix_free(gauss_A);
                    matrix_free(gauss_B);
                    matrix_free(gauss_X);
                }
                break;
                
            case 0:
                printf("Выход...\n");
                break;
                
            default:
                printf("Неверный выбор!\n");
        }
    } while(choice != 0);

    matrix_free(A);
    matrix_free(B);
    
    return 0;
}