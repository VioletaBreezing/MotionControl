#ifndef __LIE_SOPHUS_H__
#define __LIE_SOPHUS_H__

#include "MatrixMath.h"

#ifdef MATRIXMATH_USE_FLOAT
#define LIE_SOPHUS_USE_FLOAT MATRIXMATH_USE_FLOAT
#endif

// 精度类型定义
#if defined(LIE_SOPHUS_USE_FLOAT) && LIE_SOPHUS_USE_FLOAT == 1
typedef float lie_scalar_t;
#else
typedef double lie_scalar_t;
#endif

#ifdef __cplusplus
extern "C" {
#endif

// Lie代数操作函数
void skew(const lie_scalar_t* phi, lie_scalar_t* S);
void SO3_to_so3(const lie_scalar_t* R, lie_scalar_t* phi);
void so3_to_SO3(const lie_scalar_t* phi, lie_scalar_t* R);
void J_left(const lie_scalar_t* phi, lie_scalar_t* J);
void SE3_to_se3(const lie_scalar_t* T, lie_scalar_t* xi);
void se3_to_SE3(const lie_scalar_t* xi, lie_scalar_t* T);

// 辅助函数
lie_scalar_t vector_norm(const lie_scalar_t* v, int n);
void vector_normalize(const lie_scalar_t* v, int n, lie_scalar_t* result);
void matrix_add_scaled(lie_scalar_t* A, lie_scalar_t* B, lie_scalar_t alpha, int m, int n, lie_scalar_t* C);
void outer_product(const lie_scalar_t* a, const lie_scalar_t* b, int m, int n, lie_scalar_t* result);

#ifdef __cplusplus
}
#endif

#endif
