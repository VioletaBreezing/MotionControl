#ifndef __OPTIMIZER_LM_H__
#define __OPTIMIZER_LM_H__

#include "MatrixMath.h"

// 精度类型定义（与LieSophus保持一致）
#if defined(MATRIXMATH_USE_FLOAT) && MATRIXMATH_USE_FLOAT == 1
typedef float opt_scalar_t;
#else
typedef double opt_scalar_t;
#endif

// 内存分配控制宏
// 如果定义了OPT_LM_FORBIDDEN_DYNAMIC_ALLOCATION且为1，则禁用动态内存分配，使用静态内存
#ifndef OPT_LM_FORBIDDEN_DYNAMIC_ALLOCATION
#define OPT_LM_FORBIDDEN_DYNAMIC_ALLOCATION 0
#endif

// 静态内存模式下的最大尺寸定义
#if OPT_LM_FORBIDDEN_DYNAMIC_ALLOCATION
#ifndef OPT_LM_MAX_PARAM_SIZE
#define OPT_LM_MAX_PARAM_SIZE 20
#endif
#ifndef OPT_LM_MAX_RESIDUAL_SIZE
#define OPT_LM_MAX_RESIDUAL_SIZE 100
#endif
#endif

#ifdef __cplusplus
extern "C" {
#endif

// 函数指针类型定义
// 残差函数: 输入参数x和额外参数args，输出残差向量r
typedef void (*residual_func_t)(const opt_scalar_t* x, void* args, opt_scalar_t* r);

// 雅可比矩阵函数: 输入参数x和额外参数args，输出雅可比矩阵J
typedef void (*jacobian_func_t)(const opt_scalar_t* x, void* args, opt_scalar_t* J);

// LM优化器结果结构
typedef struct {
    opt_scalar_t* x_hat;        // 优化后的参数
    opt_scalar_t* resnorm;      // 最终残差向量
    int converged;              // 是否收敛 (0: false, 1: true)
    int iter;                   // 实际迭代次数
    opt_scalar_t final_cost;     // 最终代价函数值
} lm_result_t;

// LM优化器主函数
int Optimizer_LM(
    residual_func_t ResFcn,     // 残差函数指针
    void* Rargs,                // 残差函数额外参数
    jacobian_func_t JacobiFcn,  // 雅可比函数指针  
    void* Jargs,                // 雅可比函数额外参数
    const opt_scalar_t* x0,     // 初始参数估计
    int param_size,             // 参数维度 (对应MATLAB中的6)
    int residual_size,          // 残差维度
    opt_scalar_t TolX,          // 参数变化容差
    opt_scalar_t TolR,          // 残差容差  
    opt_scalar_t TolLamda,      // lambda上限容差
    int MaxIter,                // 最大迭代次数
    int debug,                  // 调试标志 (0: 关闭, 1: 开启)
    lm_result_t* result         // 输出结果结构
);

// 辅助函数
opt_scalar_t vector_norm(const opt_scalar_t* v, int n);
opt_scalar_t vector_dot(const opt_scalar_t* a, const opt_scalar_t* b, int n);

#ifdef __cplusplus
}
#endif

#endif