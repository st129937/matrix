#include "matrix.h"
#include "matrix_ops.h"
#include "matrix_manip.h"
#include "matrix_algorithms.h"
#include <stdio.h>
#include <locale.h>

void print_menu() {
    printf("\n=== ���� �������� � ��������� ===\n");
    printf("1. �������� �������� �������\n");
    printf("2. ������� �������\n");
    printf("3. ������� �������\n");
    printf("4. �������� �� ������\n");
    printf("5. ��������������� �������\n");
    printf("6. �������� ������ �������\n");
    printf("7. �������� ������� �������\n");
    printf("8. �������� �������\n");
    printf("9. ����������� �������\n");
    printf("10. ��������� ����������\n");
    printf("0. �����\n");
    printf("�������� ��������: ");
}

void print_matrix_choice() {
    printf("\n�������� �������:\n");
    printf("1. ������� A\n");
    printf("2. ������� B\n");
    printf("��� �����: ");
}

int main() {
    setlocale(LC_CTYPE, "Russian");
    
    matrix *A = matrix_alloc(2, 2);
    *matrix_ptr(A, 0, 0) = 1; *matrix_ptr(A, 0, 1) = 2;
    *matrix_ptr(A, 1, 0) = 3; *matrix_ptr(A, 1, 1) = 4;

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
                printf("\n������� A:\n");
                matrix_print(A);
                printf("\n������� B:\n");
                matrix_print(B);
                break;
                
            case 2: 
                result = matrix_copy(A);
                if(matrix_add(result, B)) {
                    printf("������: ������ ������� ������!\n");
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
                    printf("������: ������ ������� ������!\n");
                    matrix_free(result);
                } else {
                    printf("\nA - B =\n");
                    matrix_print(result);
                    matrix_free(result);
                }
                break;
                
            case 4: 
                printf("������� ������: ");
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
                printf("\n����� ������������ �����:\n");
                matrix_print(result);
                matrix_free(result);
                break;
                
            case 7: 
                result = matrix_copy(A);
                matrix_swap_cols(result, 0, 1);
                printf("\n����� ������������ ��������:\n");
                matrix_print(result);
                matrix_free(result);
                break;
                
            case 8: 
                result = matrix_alloc(2, 2);
                if(matrix_mul(result, A, B)) {
                    printf("������ ���������!\n");
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
                    printf("\n����� ������� A:\n");
                } else if(m_choice == 2) {
                    result = matrix_copy(B);
                    printf("\n����� ������� B:\n");
                } else {
                    result = NULL;
                    printf("�������� �����!\n");
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
                    printf("������ ���������� ����������!\n");
                }
                break;
                
            case 0:
                printf("�����...\n");
                break;
                
            default:
                printf("�������� �����!\n");
        }
    } while(choice != 0);

    matrix_free(A);
    matrix_free(B);
    
    return 0;
}