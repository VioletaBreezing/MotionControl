#ifndef MatrixMath_h
#define MatrixMath_h

#ifdef SOC_6678
#include "ti/dsplib/dsplib.h"
#endif

#if defined(MATRIXMATH_USE_FLOAT) && MATRIXMATH_USE_FLOAT == 1
typedef float mtx_type;
#define mtx_sin sinf
#define mtx_cos cosf
#define mtx_sqrt sqrtf
#define mtx_fabs fabsf
#else
typedef double mtx_type;
#define mtx_sin sin
#define mtx_cos cos
#define mtx_sqrt sqrt
#define mtx_fabs fabs
#endif

#ifdef __cplusplus
extern "C" {
#endif

// Function declarations for matrix operations
void MatrixMath_Print(mtx_type* A, int m, int n, const char* label);
void MatrixMath_Copy(mtx_type* A, int n, int m, mtx_type* B);
void MatrixMath_Multiply(mtx_type* A, mtx_type* B, int m, int p, int n, mtx_type* C);
void MatrixMath_Add(mtx_type* A, mtx_type* B, int m, int n, mtx_type* C);
void MatrixMath_Subtract(mtx_type* A, mtx_type* B, int m, int n, mtx_type* C);
void MatrixMath_Transpose(mtx_type* A, int m, int n, mtx_type* C);
void MatrixMath_Scale(mtx_type* A, int m, int n, mtx_type k);
int MatrixMath_Invert(mtx_type* A, int n);

#ifdef __cplusplus
}
#endif

#endif
