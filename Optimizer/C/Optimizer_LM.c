#include "Optimizer_LM.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// 数学函数适配（根据精度宏自动选择）
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

// 静态内存缓冲区（仅在禁用动态分配时使用）
#if OPT_LM_FORBIDDEN_DYNAMIC_ALLOCATION
static opt_scalar_t x_buffer[OPT_LM_MAX_PARAM_SIZE];
static opt_scalar_t r_buffer[OPT_LM_MAX_RESIDUAL_SIZE];
static opt_scalar_t J_buffer[OPT_LM_MAX_RESIDUAL_SIZE * OPT_LM_MAX_PARAM_SIZE];
static opt_scalar_t dx_buffer[OPT_LM_MAX_PARAM_SIZE];
static opt_scalar_t x_new_buffer[OPT_LM_MAX_PARAM_SIZE];
static opt_scalar_t r_new_buffer[OPT_LM_MAX_RESIDUAL_SIZE];
static opt_scalar_t JT_buffer[OPT_LM_MAX_PARAM_SIZE * OPT_LM_MAX_RESIDUAL_SIZE];
static opt_scalar_t JTJ_buffer[OPT_LM_MAX_PARAM_SIZE * OPT_LM_MAX_PARAM_SIZE];
static opt_scalar_t JTr_buffer[OPT_LM_MAX_PARAM_SIZE];
static opt_scalar_t A_buffer[OPT_LM_MAX_PARAM_SIZE * OPT_LM_MAX_PARAM_SIZE];
static opt_scalar_t A_copy_buffer[OPT_LM_MAX_PARAM_SIZE * OPT_LM_MAX_PARAM_SIZE];
static opt_scalar_t dx_new_buffer[OPT_LM_MAX_PARAM_SIZE];
#endif

// 辅助函数实现（向量操作保留，因为MatrixMath没有提供）
opt_scalar_t vector_norm(const opt_scalar_t* v, int n) {
    opt_scalar_t sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += v[i] * v[i];
    }
    return opt_sqrt(sum);
}

opt_scalar_t vector_dot(const opt_scalar_t* a, const opt_scalar_t* b, int n) {
    opt_scalar_t sum = 0.0;
    for (int i = 0; i < n; i++) {
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

    printf("[%s]: LM Optimizer Using Static Memory.\n", __func__);
    
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
    opt_scalar_t lambda = 1e-3;
    opt_scalar_t nu = 2.0;
    result->converged = 0;
    result->iter = 0;
    result->final_cost = 0.0;
    
    // 主迭代循环
    for (int iter = 1; iter <= MaxIter; iter++) {
        // 计算残差和雅可比矩阵
        ResFcn(x, Rargs, r);
        JacobiFcn(x, Jargs, J);
        
        // 计算代价函数
        opt_scalar_t cost = 0.5 * vector_dot(r, r, residual_size);
        
        // 构建正规方程: (J'J + lambda*I) * dx = -J'r
        // 使用MatrixMath库进行矩阵运算
        
#if OPT_LM_FORBIDDEN_DYNAMIC_ALLOCATION
        opt_scalar_t* JT = JT_buffer;
        opt_scalar_t* JTJ = JTJ_buffer;
        opt_scalar_t* JTr = JTr_buffer;
        opt_scalar_t* A = A_buffer;
        opt_scalar_t* A_copy = A_copy_buffer;
        opt_scalar_t* dx_new = dx_new_buffer;
#else
        opt_scalar_t* JT = (opt_scalar_t*)malloc(param_size * residual_size * sizeof(opt_scalar_t));
        opt_scalar_t* JTJ = (opt_scalar_t*)malloc(param_size * param_size * sizeof(opt_scalar_t));
        opt_scalar_t* JTr = (opt_scalar_t*)malloc(param_size * sizeof(opt_scalar_t));
        opt_scalar_t* A = (opt_scalar_t*)malloc(param_size * param_size * sizeof(opt_scalar_t));
        opt_scalar_t* A_copy = (opt_scalar_t*)malloc(param_size * param_size * sizeof(opt_scalar_t));
        opt_scalar_t* dx_new = (opt_scalar_t*)malloc(param_size * sizeof(opt_scalar_t));
        
        if (!JT || !JTJ || !JTr || !A || !A_copy || !dx_new) {
            FREE_PTR(JT); FREE_PTR(JTJ); FREE_PTR(JTr); FREE_PTR(A); FREE_PTR(A_copy); FREE_PTR(dx_new);
            FREE_PTR(x); FREE_PTR(r); FREE_PTR(J); FREE_PTR(dx); FREE_PTR(x_new); FREE_PTR(r_new);
            FREE_PTR(result->x_hat); FREE_PTR(result->resnorm);
            return -2;
        }
#endif
        
        // 计算 J' (转置)
        MatrixMath_Transpose(J, residual_size, param_size, JT);
        
        // 计算 JTJ = J' * J
        MatrixMath_Multiply(JT, J, param_size, residual_size, param_size, JTJ);
        
        // 计算 JTr = J' * r
        MatrixMath_Multiply(JT, r, param_size, residual_size, 1, JTr);
        
        // A = JTJ + lambda * I
        memcpy(A, JTJ, param_size * param_size * sizeof(opt_scalar_t));
        // 添加lambda到对角线
        for (int i = 0; i < param_size; i++) {
            A[i * param_size + i] += lambda;
        }
        
        // b = -JTr
        for (int i = 0; i < param_size; i++) {
            dx[i] = -JTr[i];
        }
        
        // 求解线性方程组 A * dx = b
        memcpy(A_copy, A, param_size * param_size * sizeof(opt_scalar_t));
        
        // 使用MatrixMath的矩阵求逆
        if (MatrixMath_Invert(A_copy, param_size) != 0) {
            // 矩阵不可逆
            if (debug) {
                printf("[%s]: Matrix inversion failed at iteration %d\n", __func__, iter);
            }
#if !OPT_LM_FORBIDDEN_DYNAMIC_ALLOCATION
            FREE_PTR(dx_new); FREE_PTR(A_copy); FREE_PTR(A); FREE_PTR(JT); FREE_PTR(JTJ); FREE_PTR(JTr);
            FREE_PTR(x); FREE_PTR(r); FREE_PTR(J); FREE_PTR(dx); FREE_PTR(x_new); FREE_PTR(r_new);
            FREE_PTR(result->x_hat); FREE_PTR(result->resnorm);
#endif
            return -3;
        }
        
        // dx = A^(-1) * b
        MatrixMath_Multiply(A_copy, dx, param_size, param_size, 1, dx_new);
        memcpy(dx, dx_new, param_size * sizeof(opt_scalar_t));
        
        // 调试输出
        if (debug) {
            printf("[%s]: Iter %2d: cost = %.3e\n", __func__, iter, cost);
        }
        
        // 收敛判断1: 残差足够小
        if (vector_norm(r, residual_size) < TolR) {
            if (debug) {
                printf("[%s]: Converged by residual norm.\n", __func__);
            }
            result->converged = 1;
            result->iter = iter;
            result->final_cost = cost;
            memcpy(result->x_hat, x, param_size * sizeof(opt_scalar_t));
            memcpy(result->resnorm, r, residual_size * sizeof(opt_scalar_t));
            
            // 清理内存
#if !OPT_LM_FORBIDDEN_DYNAMIC_ALLOCATION
            FREE_PTR(dx_new); FREE_PTR(A_copy); FREE_PTR(A); FREE_PTR(JT); FREE_PTR(JTJ); FREE_PTR(JTr);
            FREE_PTR(x); FREE_PTR(r); FREE_PTR(J); FREE_PTR(dx); FREE_PTR(x_new); FREE_PTR(r_new);
#endif
            return 0; // 成功收敛
        }
        
        // 收敛判断2: 参数变化足够小
        if (vector_norm(dx, param_size) < TolX) {
            if (vector_norm(r, residual_size) > TolR) {
                if (debug) {
                    printf("[%s]: LM failed: Residual(%.3e) norm cannot converge.\n", __func__, vector_norm(r, residual_size));
                }
            }
            if (debug) {
                printf("[%s]: Ended iteration by X change.\n", __func__);
            }
            result->converged = 0;
            result->iter = iter;
            result->final_cost = cost;
            memcpy(result->x_hat, x, param_size * sizeof(opt_scalar_t));
            memcpy(result->resnorm, r, residual_size * sizeof(opt_scalar_t));
            
            // 清理内存
#if !OPT_LM_FORBIDDEN_DYNAMIC_ALLOCATION
            FREE_PTR(dx_new); FREE_PTR(A_copy); FREE_PTR(A); FREE_PTR(JT); FREE_PTR(JTJ); FREE_PTR(JTr);
            FREE_PTR(x); FREE_PTR(r); FREE_PTR(J); FREE_PTR(dx); FREE_PTR(x_new); FREE_PTR(r_new);
#endif
            return 0; // 正常结束但未完全收敛
        }
        
        // 尝试新参数
        for (int i = 0; i < param_size; i++) {
            x_new[i] = x[i] + dx[i];
        }
        
        ResFcn(x_new, Rargs, r_new);
        opt_scalar_t cost_new = 0.5 * vector_dot(r_new, r_new, residual_size);
        
        // 计算rho
        opt_scalar_t numerator = cost - cost_new;
        opt_scalar_t denominator = 0.5 * vector_dot(dx, dx, param_size) * lambda - 
                                  0.5 * vector_dot(dx, JTr, param_size);
        opt_scalar_t rho = (denominator > OPT_EPS) ? (numerator / denominator) : 0.0;
        
        // 接受或拒绝更新
        if (rho > 0) {
            // 接受更新
            memcpy(x, x_new, param_size * sizeof(opt_scalar_t));
            lambda = lambda * ((2 * rho - 1) > 0 ? 
                              opt_fmax(1.0/3.0, 1.0 - opt_pow(2 * rho - 1, 3)) : 1.0/3.0);
            nu = 2.0;
        } else {
            // 拒绝更新
            lambda = lambda * nu;
            nu = 2.0 * nu;
        }
        
        // lambda过大检查
        if (lambda > TolLamda) {
            if (debug) {
                printf("[%s]: failed: lambda = %.4e too large.\n", __func__, lambda);
            }
            result->converged = 0;
            result->iter = iter;
            result->final_cost = cost;
            memcpy(result->x_hat, x, param_size * sizeof(opt_scalar_t));
            memcpy(result->resnorm, r, residual_size * sizeof(opt_scalar_t));
            
            // 清理内存
#if !OPT_LM_FORBIDDEN_DYNAMIC_ALLOCATION
            FREE_PTR(dx_new); FREE_PTR(A_copy); FREE_PTR(A); FREE_PTR(JT); FREE_PTR(JTJ); FREE_PTR(JTr);
            FREE_PTR(x); FREE_PTR(r); FREE_PTR(J); FREE_PTR(dx); FREE_PTR(x_new); FREE_PTR(r_new);
#endif
            return -4; // lambda过大失败
        }
        
        // 清理临时内存
#if !OPT_LM_FORBIDDEN_DYNAMIC_ALLOCATION
        FREE_PTR(dx_new); FREE_PTR(A_copy); FREE_PTR(A); FREE_PTR(JT); FREE_PTR(JTJ); FREE_PTR(JTr);
#endif
    }
    
    // 达到最大迭代次数
    result->converged = 0;
    result->iter = MaxIter;
    result->final_cost = 0.5 * vector_dot(r, r, residual_size);
    memcpy(result->x_hat, x, param_size * sizeof(opt_scalar_t));
    memcpy(result->resnorm, r, residual_size * sizeof(opt_scalar_t));
    
    // 清理内存
#if !OPT_LM_FORBIDDEN_DYNAMIC_ALLOCATION
    FREE_PTR(x); FREE_PTR(r); FREE_PTR(J); FREE_PTR(dx); FREE_PTR(x_new); FREE_PTR(r_new);
#endif
    return 0; // 达到最大迭代次数
}