#include "MatrixMath.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#if defined(MATRIXMATH_USE_FLOAT) && MATRIXMATH_USE_FLOAT == 1
#define mtx_sin sinf
#define mtx_cos cosf
#define mtx_sqrt sqrtf
#define mtx_fabs fabsf
#else
#define mtx_sin sin
#define mtx_cos cos
#define mtx_sqrt sqrt
#define mtx_fabs fabs
#endif

#if defined(FORBIDDEN_DYNAMIC_ALLOCATION) && FORBIDDEN_DYNAMIC_ALLOCATION == 1
#define MAX_MATRIX_SIZE (10)
static mtx_type buffer[MAX_MATRIX_SIZE * MAX_MATRIX_SIZE];
#endif

#ifdef SOC_C6678
#pragma CODE_SECTION(MatrixMath_Copy, ".sa_code")
#pragma CODE_SECTION(MatrixMath_Multiply, ".sa_code")
#pragma CODE_SECTION(MatrixMath_Add, ".sa_code")
#pragma CODE_SECTION(MatrixMath_Subtract, ".sa_code")
#pragma CODE_SECTION(MatrixMath_Transpose, ".sa_code")
#pragma CODE_SECTION(MatrixMath_Scale, ".sa_code")
#pragma CODE_SECTION(MatrixMath_Invert, ".sa_code")
#endif

void MatrixMath_Print(mtx_type* A, int m, int n, const char* label)
{
    int i, j;
    printf("%s:\n", label);
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            printf("%12.6f ", A[i * n + j]);
        }
        printf("\n");
    }
}

void MatrixMath_Copy(mtx_type* A, int m, int n, mtx_type* B)
{
    int i, j;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            B[i * n + j] = A[i * n + j];
        }
    }
}

void MatrixMath_Multiply(mtx_type* A, mtx_type* B, int m, int p, int n, mtx_type* C)
{
    int i, j, k;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            C[i * n + j] = 0;
            for (k = 0; k < p; k++) {
                C[i * n + j] += A[i * p + k] * B[k * n + j];
            }
        }
    }
}

void MatrixMath_Add(mtx_type* A, mtx_type* B, int m, int n, mtx_type* C)
{
    int i, j;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            C[i * n + j] = A[i * n + j] + B[i * n + j];
        }
    }
}

void MatrixMath_Subtract(mtx_type* A, mtx_type* B, int m, int n, mtx_type* C)
{
    int i, j;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            C[i * n + j] = A[i * n + j] - B[i * n + j];
        }
    }
}

void MatrixMath_Transpose(mtx_type* A, int m, int n, mtx_type* C)
{
    int i, j;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            C[j * m + i] = A[i * n + j];
        }
    }
}

void MatrixMath_Scale(mtx_type* A, int m, int n, mtx_type k)
{
    int i, j;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            A[i * n + j] *= k;
        }
    }
}

// 矩阵求逆函数（使用高斯-约旦消元法）
int MatrixMath_Invert(mtx_type* A, int n)
{
    int i, j, k;
    mtx_type temp;
#if defined(FORBIDDEN_DYNAMIC_ALLOCATION) && FORBIDDEN_DYNAMIC_ALLOCATION == 1
    mtx_type* I = buffer;
    if (n * n > MAX_MATRIX_SIZE * MAX_MATRIX_SIZE) {
        // 超出预定义的最大矩阵大小，无法处理
        fprintf(stderr, "Error: Matrix size exceeds maximum allowed size.\n");
        return -1;
    }   
#else
    mtx_type* I = (mtx_type*)malloc(n * n * sizeof(mtx_type));
#endif
    
    // 初始化单位矩阵
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            I[i * n + j] = (i == j) ? 1.0 : 0.0;
        }
    }
    
    // 高斯-约旦消元
    for (i = 0; i < n; i++) {
        // 寻找主元
        int pivot = i;
        mtx_type max_val = mtx_fabs(A[i * n + i]);
        for (k = i + 1; k < n; k++) {
            if (mtx_fabs(A[k * n + i]) > max_val) {
                max_val = mtx_fabs(A[k * n + i]);
                pivot = k;
            }
        }
        
        // 如果主元为0，矩阵不可逆
        if (max_val < 1e-10) {
#if defined(FORBIDDEN_DYNAMIC_ALLOCATION) && FORBIDDEN_DYNAMIC_ALLOCATION == 1
            // Do nothing, as we can't free the buffer
#else
            free(I);
#endif
            return -1;
        }
        
        // 交换行
        if (pivot != i) {
            for (j = 0; j < n; j++) {
                temp = A[i * n + j];
                A[i * n + j] = A[pivot * n + j];
                A[pivot * n + j] = temp;
                
                temp = I[i * n + j];
                I[i * n + j] = I[pivot * n + j];
                I[pivot * n + j] = temp;
            }
        }
        
        // 归一化主行
        temp = A[i * n + i];
        for (j = 0; j < n; j++) {
            A[i * n + j] /= temp;
            I[i * n + j] /= temp;
        }
        
        // 消元
        for (k = 0; k < n; k++) {
            if (k != i) {
                temp = A[k * n + i];
                for (j = 0; j < n; j++) {
                    A[k * n + j] -= temp * A[i * n + j];
                    I[k * n + j] -= temp * I[i * n + j];
                }
            }
        }
    }
    
    // 将结果复制回A
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            A[i * n + j] = I[i * n + j];
        }
    }
    
#if defined(FORBIDDEN_DYNAMIC_ALLOCATION) && FORBIDDEN_DYNAMIC_ALLOCATION == 1
    // Do nothing, as we can't free the buffer
#else
    free(I);
#endif
    return 0;
}
