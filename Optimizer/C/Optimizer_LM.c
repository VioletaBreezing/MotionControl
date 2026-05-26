#include "Optimizer_LM.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// #ifdef SOC_C6678
// #include <ti/dsplib/dsplib.h>
// #ifdef LIE_SOPHUS_USE_FLOAT

// #else
// #endif

// #else
#include <math.h>
#ifdef LIE_SOPHUS_USE_FLOAT
#define opt_sin sinf
#define opt_cos cosf
#define opt_sqrt sqrtf
#define opt_fabs fabsf
#define opt_pow powf
#define opt_fmax fmaxf
#define OPT_EPS 1e-6f
#else
#define opt_sin sin
#define opt_cos cos
#define opt_sqrt sqrt
#define opt_fabs fabs
#define opt_pow pow
#define opt_fmax fmax
#define OPT_EPS 1e-12
#endif
// #endif

#ifdef SOC_C6678
#include <ti/dsplib/dsplib.h>
#endif

// 静态内存缓冲区（仅在禁用动态分配时使用）
#if OPT_LM_FORBIDDEN_DYNAMIC_ALLOCATION
#pragma DATA_ALIGN(x_buffer, 8)
static opt_scalar_t x_buffer[OPT_LM_MAX_PARAM_SIZE];

#pragma DATA_ALIGN(r_buffer, 8)
static opt_scalar_t r_buffer[OPT_LM_MAX_RESIDUAL_SIZE];

#pragma DATA_ALIGN(J_buffer, 8)
static opt_scalar_t J_buffer[OPT_LM_MAX_RESIDUAL_SIZE * OPT_LM_MAX_PARAM_SIZE];

#pragma DATA_ALIGN(dx_buffer, 8)
static opt_scalar_t dx_buffer[OPT_LM_MAX_PARAM_SIZE];

#pragma DATA_ALIGN(x_new_buffer, 8)
static opt_scalar_t x_new_buffer[OPT_LM_MAX_PARAM_SIZE];

#pragma DATA_ALIGN(r_new_buffer, 8)
static opt_scalar_t r_new_buffer[OPT_LM_MAX_RESIDUAL_SIZE];

#pragma DATA_ALIGN(JT_buffer, 8)
static opt_scalar_t JT_buffer[OPT_LM_MAX_PARAM_SIZE * OPT_LM_MAX_RESIDUAL_SIZE];

#pragma DATA_ALIGN(L_buffer, 8)
static opt_scalar_t L_buffer[OPT_LM_MAX_PARAM_SIZE * OPT_LM_MAX_PARAM_SIZE];

#pragma DATA_ALIGN(JTr_buffer, 8)
static opt_scalar_t JTr_buffer[OPT_LM_MAX_PARAM_SIZE];

#pragma DATA_ALIGN(A_buffer, 8)
static opt_scalar_t A_buffer[OPT_LM_MAX_PARAM_SIZE * OPT_LM_MAX_PARAM_SIZE];

#pragma DATA_ALIGN(A_copy_buffer, 8)
static opt_scalar_t A_copy_buffer[OPT_LM_MAX_PARAM_SIZE * OPT_LM_MAX_PARAM_SIZE];

#pragma DATA_ALIGN(dx_new_buffer, 8)
static opt_scalar_t dx_new_buffer[OPT_LM_MAX_PARAM_SIZE];

#pragma DATA_ALIGN(b_buffer, 8)
static opt_scalar_t b_buffer [OPT_LM_MAX_PARAM_SIZE];

#pragma DATA_ALIGN(y_buffer, 8)
static opt_scalar_t y_buffer [OPT_LM_MAX_PARAM_SIZE];
#endif

#pragma DATA_ALIGN(dx_new_buffer, 8)
static opt_scalar_t EYE6[36] = {1,0,0,0,0,0,
                                0,1,0,0,0,0,
                                0,0,1,0,0,0,
                                0,0,0,1,0,0,
                                0,0,0,0,1,0,
                                0,0,0,0,0,1};

#ifdef SOC_C6678
#pragma CODE_SECTION(opt_vector_norm, ".sa_code")
#pragma CODE_SECTION(opt_vector_dot, ".sa_code")
#pragma CODE_SECTION(Optimizer_LM, ".sa_code")
#endif

// 辅助函数实现（向量操作保留，因为MatrixMath没有提供）
opt_scalar_t opt_vector_norm(const opt_scalar_t* v, int n) {
    opt_scalar_t sum = 0.0;
    int i;
    for (i = 0; i < n; i++) {
        sum += v[i] * v[i];
    }
    return opt_sqrt(sum);
}

opt_scalar_t opt_vector_dot(const opt_scalar_t* a, const opt_scalar_t* b, int n) {
    opt_scalar_t sum = 0.0;
    int i;
    for (i = 0; i < n; i++) {
        sum += a[i] * b[i];
    }
    return sum;
}

// 内存分配辅助宏
#if OPT_LM_FORBIDDEN_DYNAMIC_ALLOCATION
#define ALLOC_VECTOR(size, max_size, buffer) \
    ((size) <= (max_size) ? (buffer) : NULL)
#define ALLOC_MATRIX(rows, cols, max_rows, max_cols, buffer) \
    (((rows) <= (max_rows) && (cols) <= (max_cols)) ? (buffer) : NULL)
#define FREE_PTR(ptr) /* 静态模式下不需要释放 */
#else
#define ALLOC_VECTOR(size, max_size, buffer) \
    ((opt_scalar_t*)malloc((size) * sizeof(opt_scalar_t)))
#define ALLOC_MATRIX(rows, cols, max_rows, max_cols, buffer) \
    ((opt_scalar_t*)malloc((rows) * (cols) * sizeof(opt_scalar_t)))
#define FREE_PTR(ptr) free(ptr)
#endif

// LM优化器主函数
int Optimizer_LM(
    residual_func_t ResFcn,
    void* Rargs,
    jacobian_func_t JacobiFcn,
    void* Jargs,
    const opt_scalar_t* x0,
    int param_size,
    int residual_size,
    opt_scalar_t TolX,
    opt_scalar_t TolR,
    opt_scalar_t TolLamda,
    int MaxIter,
    int debug,
    lm_result_t* result) {
    
    // 输入参数验证
    if (!ResFcn || !JacobiFcn || !x0 || !result || param_size <= 0 || residual_size <= 0) {
        return -1; // 参数错误
    }
    
#if OPT_LM_FORBIDDEN_DYNAMIC_ALLOCATION
    // 静态内存分配模式：检查尺寸限制
    if (param_size > OPT_LM_MAX_PARAM_SIZE || residual_size > OPT_LM_MAX_RESIDUAL_SIZE) {
        return -5; // 尺寸超出静态缓冲区限制
    }
    
    // 使用静态缓冲区
    opt_scalar_t* x = x_buffer;
    opt_scalar_t* r = r_buffer;
    opt_scalar_t* J = J_buffer;
    opt_scalar_t* dx = dx_buffer;
    opt_scalar_t* x_new = x_new_buffer;
    opt_scalar_t* r_new = r_new_buffer;
    
    // 在静态模式下，结果内存必须由调用者提供
    if (!result->x_hat || !result->resnorm) {
        return -1; // 结果内存未分配
    }

    // printf("[%s]: LM Optimizer Using Static Memory.\n", __func__);
    
#else
    // 动态内存分配模式
    opt_scalar_t* x = (opt_scalar_t*)malloc(param_size * sizeof(opt_scalar_t));
    opt_scalar_t* r = (opt_scalar_t*)malloc(residual_size * sizeof(opt_scalar_t));
    opt_scalar_t* J = (opt_scalar_t*)malloc(residual_size * param_size * sizeof(opt_scalar_t));
    opt_scalar_t* dx = (opt_scalar_t*)malloc(param_size * sizeof(opt_scalar_t));
    opt_scalar_t* x_new = (opt_scalar_t*)malloc(param_size * sizeof(opt_scalar_t));
    opt_scalar_t* r_new = (opt_scalar_t*)malloc(residual_size * sizeof(opt_scalar_t));
    
    // 结果内存分配
    result->x_hat = (opt_scalar_t*)malloc(param_size * sizeof(opt_scalar_t));
    result->resnorm = (opt_scalar_t*)malloc(residual_size * sizeof(opt_scalar_t));
    
    if (!x || !r || !J || !dx || !x_new || !r_new || !result->x_hat || !result->resnorm) {
        // 内存分配失败，清理并返回错误
        FREE_PTR(x); FREE_PTR(r); FREE_PTR(J); FREE_PTR(dx); FREE_PTR(x_new); FREE_PTR(r_new);
        if (result->x_hat) FREE_PTR(result->x_hat);
        if (result->resnorm) FREE_PTR(result->resnorm);
        return -2; // 内存分配失败
    }
#endif
    
    // 初始化
    memcpy(x, x0, param_size * sizeof(opt_scalar_t));
    opt_scalar_t lambda = 1e-5;
    opt_scalar_t nu = 2.0;
    result->converged = 0;
    result->iter = 0;
    result->final_cost = 0.0;

    // 主迭代循环
    int iter;
    for (iter = 1; iter <= MaxIter; iter++) {
        int i;

        // 计算残差和雅可比矩阵
        ResFcn(x, Rargs, r);
        JacobiFcn(x, Jargs, J);
        
        // 计算代价函数
        opt_scalar_t cost = 0.5 * opt_vector_dot(r, r, residual_size);
        
        // 构建正规方程: (J'J + lambda*I) * dx = -J'r
        // 使用MatrixMath库进行矩阵运算
        
#if OPT_LM_FORBIDDEN_DYNAMIC_ALLOCATION
        opt_scalar_t* JT = JT_buffer;
        opt_scalar_t* L = L_buffer;
        opt_scalar_t* JTr = JTr_buffer;
        opt_scalar_t* A = A_buffer;
        opt_scalar_t* b = b_buffer;
        opt_scalar_t* y = y_buffer;
        opt_scalar_t* A_copy = A_copy_buffer;
        opt_scalar_t* dx_new = dx_new_buffer;
#else
        opt_scalar_t* JT = (opt_scalar_t*)malloc(param_size * residual_size * sizeof(opt_scalar_t));
        opt_scalar_t* L = (opt_scalar_t*)malloc(param_size * param_size * sizeof(opt_scalar_t));
        opt_scalar_t* JTr = (opt_scalar_t*)malloc(param_size * sizeof(opt_scalar_t));
        opt_scalar_t* A = (opt_scalar_t*)malloc(param_size * param_size * sizeof(opt_scalar_t));
        opt_scalar_t* b = (opt_scalar_t*)malloc(param_size * sizeof(opt_scalar_t));
        opt_scalar_t* y = (opt_scalar_t*)malloc(param_size * sizeof(opt_scalar_t));
        opt_scalar_t* A_copy = (opt_scalar_t*)malloc(param_size * param_size * sizeof(opt_scalar_t));
        opt_scalar_t* dx_new = (opt_scalar_t*)malloc(param_size * sizeof(opt_scalar_t));
        
        if (!JT || !L || !JTr || !A || !A_copy || !dx_new || !b || !y) {
            FREE_PTR(JT); FREE_PTR(L); FREE_PTR(JTr); FREE_PTR(A); FREE_PTR(A_copy); FREE_PTR(dx_new);
            FREE_PTR(x); FREE_PTR(r); FREE_PTR(J); FREE_PTR(dx); FREE_PTR(x_new); FREE_PTR(r_new);
            FREE_PTR(result->x_hat); FREE_PTR(result->resnorm); FREE_PTR(b); FREE_PTR(y);
            return -2;
        }
#endif
        memset(JT, 0x00, param_size * residual_size * sizeof(opt_scalar_t));
        memset(L, 0x00, param_size * param_size * sizeof(opt_scalar_t));
        memset(JTr, 0x00, param_size * sizeof(opt_scalar_t));
        memset(A, 0x00, param_size * param_size * sizeof(opt_scalar_t));
        memset(A_copy, 0x00, param_size * param_size * sizeof(opt_scalar_t));
        memset(dx_new, 0x00, param_size * sizeof(opt_scalar_t));
        
        // 计算 J' (转置)
        #ifdef SOC_6678
        DSPF_dp_mat_trans(J, residual_size, param_size, JT);
        #else
        MatrixMath_Transpose(J, residual_size, param_size, JT);
        #endif
        
        // 计算 A = JTJ = J' * J
        // MatrixMath_Multiply(JT, J, param_size, residual_size, param_size, A);
        #ifdef SOC_6678
        DSPF_dp_mat_mul_gemm_cn(JT, 1.0, param_size, residual_size, J, param_size, A);
        #else
        MatrixMath_Multiply(JT, J, param_size, residual_size, param_size, A);
        #endif
        
        // 计算 JTr = J' * r
        MatrixMath_Multiply(JT, r, param_size, residual_size, 1, JTr);
        // DSPF_dp_mat_mul_gemm_cn(JT, 1.0, param_size, residual_size, r, 1, JTr);
        
        // A = A + lambda * I
        // memcpy(A, L, param_size * param_size * sizeof(opt_scalar_t));
        // 添加lambda到对角线
        
        for (i = 0; i < param_size; i++) {
            A[i * param_size + i] += lambda;
        }
        // A[0] += lambda;  A[7] += lambda;  A[14] += lambda;
        // A[21] += lambda; A[28] += lambda; A[35] += lambda;
        
        // b = -JTr
        // for (i = 0; i < param_size; i++) {
        //     dx[i] = -JTr[i];
        // }
        b[0] = -JTr[0]; b[1] = -JTr[1]; b[2] = -JTr[2];
        b[3] = -JTr[3]; b[4] = -JTr[4]; b[5] = -JTr[5];
        
        // 求解线性方程组 A * dx = b
        memcpy(A_copy, A, param_size * param_size * sizeof(opt_scalar_t));
        
        // 使用MatrixMath的矩阵求逆
        if (MatrixMath_Invert(A_copy, param_size) != 0) {
            // 矩阵不可逆
            if (debug) {
                printf("Optimizer_LM: Matrix inversion failed at iteration %d\n", iter);
            }
#if !OPT_LM_FORBIDDEN_DYNAMIC_ALLOCATION
            FREE_PTR(dx_new); FREE_PTR(A_copy); FREE_PTR(A); FREE_PTR(JT); FREE_PTR(L); FREE_PTR(JTr);
            FREE_PTR(x); FREE_PTR(r); FREE_PTR(J); FREE_PTR(dx); FREE_PTR(x_new); FREE_PTR(r_new);
            FREE_PTR(result->x_hat); FREE_PTR(result->resnorm); FREE_PTR(b); FREE_PTR(y);
#endif
            return -3;
        }
        
        // dx = A^(-1) * b
        MatrixMath_Multiply(A_copy, b, param_size, param_size, 1, dx_new);
        memcpy(dx, dx_new, param_size * sizeof(opt_scalar_t));

        #ifdef SOC_6678
        // Cholesky 分解
        DSPF_dp_cholesky(0, param_size, A, L);
        // Cholesky 解线性方程 A*dx = b; L*LT*dx=b; L*y=b;
        DSPF_dp_cholesky_solver(param_size, L, y, b, dx);
        #endif
        
        // 调试输出
        if (debug) {
            printf("[Optimizer_LM] Iter %2d: cost = %.3e\n", iter, cost);
        }
        
        // 收敛判断2: 参数变化足够小
        if (/*opt_vector_norm(dx, param_size) < TolX && */opt_vector_norm(r, residual_size) < TolR) {
            if (opt_vector_norm(r, residual_size) > TolR) {
                if (debug) {
                    printf("[Optimizer_LM] LM failed: Residual(%.3e) norm cannot converge.\n", 
                           opt_vector_norm(r, residual_size));
                }
            }
            if (debug) {
                printf("[Optimizer_LM]: Ended iteration by X change.\n");
            }
            result->converged = 0;
            result->iter = iter;
            result->final_cost = cost;
            memcpy(result->x_hat, x, param_size * sizeof(opt_scalar_t));
            memcpy(result->resnorm, r, residual_size * sizeof(opt_scalar_t));
            
            // 清理内存
#if !OPT_LM_FORBIDDEN_DYNAMIC_ALLOCATION
            FREE_PTR(dx_new); FREE_PTR(A_copy); FREE_PTR(A); FREE_PTR(JT); FREE_PTR(L); FREE_PTR(JTr);
            FREE_PTR(x); FREE_PTR(r); FREE_PTR(J); FREE_PTR(dx); FREE_PTR(x_new); FREE_PTR(r_new);
            FREE_PTR(b); FREE_PTR(y);
#endif
            return 0; // 正常结束但未完全收敛
        }
        
        // 尝试新参数
        for (i = 0; i < param_size; i++) {
            x_new[i] = x[i] + dx[i];
        }
        
        ResFcn(x_new, Rargs, r_new);
        opt_scalar_t cost_new = 0.5 * opt_vector_dot(r_new, r_new, residual_size);
        
        // 计算rho
        opt_scalar_t numerator = cost - cost_new;
        opt_scalar_t Jdx[residual_size];
        MatrixMath_Multiply(J, dx, residual_size, param_size, 1, Jdx); // J*dx
        opt_scalar_t JTr_dx = opt_vector_dot(r, Jdx, residual_size);
        opt_scalar_t denominator = - JTr_dx - 0.5 * opt_vector_dot(Jdx, Jdx, residual_size);
        // opt_scalar_t denominator = 0.5 * (opt_vector_dot(r, r, residual_size) * lambda - opt_vector_dot(dx, JTr, param_size));
        // opt_scalar_t rho = (fabs(denominator) > OPT_EPS) ? (numerator / denominator) : 0.0;
        opt_scalar_t rho = numerator / denominator;
        
        // 接受或拒绝更新
        if (rho > 1e-4) {
            // 接受更新
            memcpy(x, x_new, param_size * sizeof(opt_scalar_t));
            lambda = lambda * ((2 * rho - 1) > 0 ? opt_fmax(1.0/3.0, 1.0 - opt_pow(2 * rho - 1, 3)) : 1.0/3.0);
            // lambda = lambda / nu;
            nu = 2.0;
        } else{
            // 拒绝更新
            lambda = lambda * nu;
            nu = 2.0 * nu;
        }
        
        // lambda过大检查
        if (lambda > TolLamda) {
            if (debug) {
                printf("[Optimizer_LM] failed: lambda = %.4e too large.\n", lambda);
            }
            result->converged = 0;
            result->iter = iter;
            result->final_cost = cost;
            memcpy(result->x_hat, x, param_size * sizeof(opt_scalar_t));
            memcpy(result->resnorm, r, residual_size * sizeof(opt_scalar_t));
            
            // 清理内存
#if !OPT_LM_FORBIDDEN_DYNAMIC_ALLOCATION
            FREE_PTR(dx_new); FREE_PTR(A_copy); FREE_PTR(A); FREE_PTR(JT); FREE_PTR(L); FREE_PTR(JTr);
            FREE_PTR(x); FREE_PTR(r); FREE_PTR(J); FREE_PTR(dx); FREE_PTR(x_new); FREE_PTR(r_new);
#endif
            return -4; // lambda过大失败
        }
        
        // 清理临时内存
#if !OPT_LM_FORBIDDEN_DYNAMIC_ALLOCATION
        FREE_PTR(dx_new); FREE_PTR(A_copy); FREE_PTR(A); FREE_PTR(JT); FREE_PTR(L); FREE_PTR(JTr);
#endif
    }
    
    // 达到最大迭代次数
    result->converged = 0;
    result->iter = MaxIter;
    result->final_cost = 0.5 * opt_vector_dot(r, r, residual_size);
    memcpy(result->x_hat, x, param_size * sizeof(opt_scalar_t));
    memcpy(result->resnorm, r, residual_size * sizeof(opt_scalar_t));
    
    // 清理内存
#if !OPT_LM_FORBIDDEN_DYNAMIC_ALLOCATION
    FREE_PTR(x); FREE_PTR(r); FREE_PTR(J); FREE_PTR(dx); FREE_PTR(x_new); FREE_PTR(r_new);
#endif
    return -5; // 达到最大迭代次数
}
