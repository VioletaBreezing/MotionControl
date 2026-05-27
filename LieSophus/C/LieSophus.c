#include "LieSophus.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef SOC_6678
#include "ti/dsplib/dsplib.h"
#include "ti/mathlib/mathlib.h"
#endif

// 如果 M_PI 仍未定义，则手动定义
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef LIE_SOPHUS_USE_FLOAT
#define LIE_SOPHUS_USE_FLOAT 0
#endif

// 数学函数适配（根据精度宏自动选择）
#if defined(LIE_SOPHUS_USE_FLOAT) && LIE_SOPHUS_USE_FLOAT == 1
#define lie_sin sinf
#define lie_cos cosf
#define lie_acos acosf
#define lie_sqrt sqrtf
#define lie_fabs fabsf
#else
#define lie_sin sin
#define lie_cos cos
#define lie_acos acos
#define lie_sqrt sqrt
#define lie_fabs fabs
#endif

static inline void matrix_eye3(lie_scalar_t* A);

#ifdef SOC_C6678
#pragma CODE_SECTION(vector_norm, ".sa_code")
#pragma CODE_SECTION(vector_normalize, ".sa_code")
#pragma CODE_SECTION(matrix_eye3, ".sa_code")
#pragma CODE_SECTION(matrix_add_scaled, ".sa_code")
#pragma CODE_SECTION(outer_product, ".sa_code")
#pragma CODE_SECTION(skew, ".sa_code")
#pragma CODE_SECTION(SO3_to_so3, ".sa_code")
#pragma CODE_SECTION(so3_to_SO3, ".sa_code")
#pragma CODE_SECTION(J_left, ".sa_code")
#pragma CODE_SECTION(SE3_to_se3, ".sa_code")
#pragma CODE_SECTION(se3_to_SE3, ".sa_code")
#endif

// 辅助函数实现
lie_scalar_t vector_norm(const lie_scalar_t* v, int n) {
    lie_scalar_t sum = 0.0;
    int i;
    for (i = 0; i < n; i++) {
        sum += v[i] * v[i];
    }
    return lie_sqrt(sum);
}

void vector_normalize(const lie_scalar_t* v, int n, lie_scalar_t* result) {
    lie_scalar_t norm = vector_norm(v, n);
    int i;
    if (norm > 1e-10) {
        for (i = 0; i < n; i++) {
            result[i] = v[i] / norm;
        }
    } else {
        for (i = 0; i < n; i++) {
            result[i] = 0.0;
        }
    }
}

static inline void matrix_eye3(lie_scalar_t* A) {
    A[0] = 1.0; A[1] = 0.0; A[2] = 0.0;
    A[3] = 0.0; A[4] = 1.0; A[5] = 0.0;
    A[6] = 0.0; A[7] = 0.0; A[8] = 1.0;
}

void matrix_add_scaled(lie_scalar_t* A, lie_scalar_t* B, lie_scalar_t alpha, int m, int n, lie_scalar_t* C) {
    int i, j;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            C[i * n + j] = A[i * n + j] + alpha * B[i * n + j];
        }
    }
}

void outer_product(const lie_scalar_t* a, const lie_scalar_t* b, int m, int n, lie_scalar_t* result) {
    int i, j;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            result[i * n + j] = a[i] * b[j];
        }
    }
}

// Lie代数操作函数实现

void skew(const lie_scalar_t* phi, lie_scalar_t* S) {
    // S = [   0   -phi(3)  phi(2);
    //       phi(3)    0   -phi(1);
    //      -phi(2)  phi(1)    0 ];
    S[0] = 0.0;        S[1] = -phi[2];    S[2] = phi[1];
    S[3] = phi[2];     S[4] = 0.0;        S[5] = -phi[0];
    S[6] = -phi[1];    S[7] = phi[0];     S[8] = 0.0;
}

void SO3_to_so3(const lie_scalar_t* R, lie_scalar_t* phi) {
    // 计算旋转角 theta
    lie_scalar_t cos_theta = (R[0] + R[4] + R[8] - 1.0) / 2.0;
    cos_theta = (cos_theta > 1.0) ? 1.0 : ((cos_theta < -1.0) ? -1.0 : cos_theta);
    lie_scalar_t theta = lie_acos(cos_theta);

    if (lie_fabs(theta) < 1e-13) {
        // 情况1: 无旋转 → phi = [0;0;0]
        phi[0] = 0.0;
        phi[1] = 0.0;
        phi[2] = 0.0;
    } else if (theta > M_PI - 1e-6) {
        // 情况2: 接近180度旋转（数值不稳定，需特殊处理）
        // 利用 R + I 的最大对角元确定主轴
        lie_scalar_t A[9];
        A[0] = (R[0] + 1.0) / 2.0; A[1] = R[1] / 2.0;        A[2] = R[2] / 2.0;
        A[3] = R[3] / 2.0;        A[4] = (R[4] + 1.0) / 2.0; A[5] = R[5] / 2.0;
        A[6] = R[6] / 2.0;        A[7] = R[7] / 2.0;        A[8] = (R[8] + 1.0) / 2.0;
        
        // 找到最大对角元的索引
        int idx = 0;
        lie_scalar_t max_val = A[0];
        if (A[4] > max_val) {
            max_val = A[4];
            idx = 1;
        }
        if (A[8] > max_val) {
            idx = 2;
        }
        
        lie_scalar_t n[3];
        n[0] = A[idx];
        n[1] = A[3 + idx];
        n[2] = A[6 + idx];
        
        lie_scalar_t norm_n = vector_norm(n, 3);
        if (norm_n > 1e-13) {
            n[0] /= norm_n;
            n[1] /= norm_n;
            n[2] /= norm_n;
        }
        
        phi[0] = M_PI * n[0];
        phi[1] = M_PI * n[1];
        phi[2] = M_PI * n[2];
    } else {
        // 情况3: 一般情况
        lie_scalar_t sin_theta = lie_sin(theta);
        lie_scalar_t nx = (R[7] - R[5]) / (2.0 * sin_theta);
        lie_scalar_t ny = (R[2] - R[6]) / (2.0 * sin_theta);
        lie_scalar_t nz = (R[3] - R[1]) / (2.0 * sin_theta);
        phi[0] = theta * nx;
        phi[1] = theta * ny;
        phi[2] = theta * nz;
    }
}

void so3_to_SO3(const lie_scalar_t* phi, lie_scalar_t* R) {
    lie_scalar_t theta = vector_norm(phi, 3);
    
    if (lie_fabs(theta) < 1e-13) {
        matrix_eye3(R);
    } else {
        lie_scalar_t k[3];
        k[0] = phi[0] / theta;
        k[1] = phi[1] / theta;
        k[2] = phi[2] / theta;
        
        lie_scalar_t K[9];
        skew(k, K);
        
        lie_scalar_t kkT[9];
        outer_product(k, k, 3, 3, kkT);
        
        lie_scalar_t cos_t = lie_cos(theta);
        lie_scalar_t sin_t = lie_sin(theta);
        
        // R = cos(theta)*eye(3) + (1-cos(theta))*(k*k') + sin(theta)*K;
        int i;
        for (i = 0; i < 9; i++) {
            R[i] = cos_t * ((i % 4 == 0) ? 1.0 : 0.0) + 
                   (1.0 - cos_t) * kkT[i] + 
                   sin_t * K[i];
        }
        // 修正对角线元素
        R[0] = cos_t + (1.0 - cos_t) * kkT[0];
        R[4] = cos_t + (1.0 - cos_t) * kkT[4];
        R[8] = cos_t + (1.0 - cos_t) * kkT[8];
    }
}

static inline void J_left(const lie_scalar_t* phi, lie_scalar_t* J) {
    lie_scalar_t theta = vector_norm(phi, 3);
    
    if (lie_fabs(theta) < 1e-13) {
        matrix_eye3(J);
    } else {
        lie_scalar_t a[3];
        a[0] = phi[0] / theta;
        a[1] = phi[1] / theta;
        a[2] = phi[2] / theta;
        
        lie_scalar_t A_outer[9];
        outer_product(a, a, 3, 3, A_outer);
        
        lie_scalar_t skew_a[9];
        skew(a, skew_a);
        
        lie_scalar_t sin_t = lie_sin(theta);
        lie_scalar_t cos_t = lie_cos(theta);
        
        lie_scalar_t coeff1 = sin_t / theta;
        lie_scalar_t coeff2 = 1.0 - sin_t / theta;
        lie_scalar_t coeff3 = (1.0 - cos_t) / theta;
        
        // J = sin(theta)/theta * eye(3) + (1 - sin(theta)/theta) * (a * a') + (1 - cos(theta))/theta * skew(a);
        // 直接计算每个元素，避免条件判断
        int i;
        for (i = 0; i < 9; i++) {
            lie_scalar_t eye_val = (i == 0 || i == 4 || i == 8) ? 1.0 : 0.0;
            J[i] = coeff1 * eye_val + coeff2 * A_outer[i] + coeff3 * skew_a[i];
        }
    }
}

void SE3_to_se3(const lie_scalar_t* T, lie_scalar_t* xi) {
    // 提取R和t
    lie_scalar_t R[9];
    lie_scalar_t t[3];
    
    R[0] = T[0]; R[1] = T[1]; R[2] = T[2];
    R[3] = T[4]; R[4] = T[5]; R[5] = T[6];
    R[6] = T[8]; R[7] = T[9]; R[8] = T[10];
    
    t[0] = T[3];
    t[1] = T[7];
    t[2] = T[11];
    
    // 计算phi
    lie_scalar_t phi[3];
    SO3_to_so3(R, phi);
    
    // 计算J
    lie_scalar_t J[9];
    J_left(phi, J);
    
    // 求解 rho = J \ t (即 J * rho = t)
    // 创建增广矩阵 [J | t] 并求解
    lie_scalar_t J_copy[9];
    MatrixMath_Copy(J, 3, 3, J_copy);
    
    if (MatrixMath_Invert(J_copy, 3) == 0) {
        // J可逆，计算 rho = J^(-1) * t
        xi[0] = J_copy[0] * t[0] + J_copy[1] * t[1] + J_copy[2] * t[2];
        xi[1] = J_copy[3] * t[0] + J_copy[4] * t[1] + J_copy[5] * t[2];
        xi[2] = J_copy[6] * t[0] + J_copy[7] * t[1] + J_copy[8] * t[2];
    } else {
        // J不可逆，使用最小二乘解或其他方法
        xi[0] = t[0];
        xi[1] = t[1];
        xi[2] = t[2];
    }
    
    // 设置phi部分
    xi[3] = phi[0];
    xi[4] = phi[1];
    xi[5] = phi[2];
}

void se3_to_SE3(const lie_scalar_t* xi, lie_scalar_t* T) {
    lie_scalar_t rho[3], phi[3];
    rho[0] = xi[0]; rho[1] = xi[1]; rho[2] = xi[2];
    phi[0] = xi[3]; phi[1] = xi[4]; phi[2] = xi[5];
    
    lie_scalar_t theta = vector_norm(phi, 3);
    lie_scalar_t R[9], t[3];
    
    if (lie_fabs(theta) < 1e-8) {
        // 小角度近似：R ≈ I, t ≈ rho
        matrix_eye3(R);
        t[0] = rho[0];
        t[1] = rho[1];
        t[2] = rho[2];
    } else {
        so3_to_SO3(phi, R);
        
        // 计算 t = J * rho
        lie_scalar_t J[9];
        J_left(phi, J);
        
        t[0] = J[0] * rho[0] + J[1] * rho[1] + J[2] * rho[2];
        t[1] = J[3] * rho[0] + J[4] * rho[1] + J[5] * rho[2];
        t[2] = J[6] * rho[0] + J[7] * rho[1] + J[8] * rho[2];
    }
    
    // 构造 SE(3) 矩阵    
    T[0] = R[0]; T[1] = R[1]; T[2] = R[2];  T[3] = t[0];
    T[4] = R[3]; T[5] = R[4]; T[6] = R[5];  T[7] = t[1];
    T[8] = R[6]; T[9] = R[7]; T[10] = R[8]; T[11] = t[2];
    T[12] = 0.0; T[13] = 0.0; T[14] = 0.0;  T[15] = 1.0;
}

void euler_to_SO3(const lie_scalar_t* euler, lie_scalar_t* R) {
    lie_scalar_t cx = lie_cos(euler[0]);
    lie_scalar_t sx = lie_sin(euler[0]);
    lie_scalar_t cy = lie_cos(euler[1]);
    lie_scalar_t sy = lie_sin(euler[1]);
    lie_scalar_t cz = lie_cos(euler[2]);
    lie_scalar_t sz = lie_sin(euler[2]);

    R[0] = cz * cy;
    R[1] = cz * sy * sx - sz * cx;
    R[2] = cz * sy * cx + sz * sx;
    R[3] = sz * cy;
    R[4] = sz * sy * sx + cz * cx;
    R[5] = sz * sy * cx - cz * sx;
    R[6] = -sy;
    R[7] = sx * cy;
    R[8] = cx * cy;
}