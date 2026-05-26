#ifndef __MIRROR_REFLECTION_SIMULATION_2RS_MYLM_H__
#define __MIRROR_REFLECTION_SIMULATION_2RS_MYLM_H__

#include "LieSophus.h"
#include "Optimizer_LM.h"

#ifdef __cplusplus
extern "C" {
#endif

// 光路参数结构
typedef struct {
    lie_scalar_t a[3];  // 入射光线起点
    lie_scalar_t u[3];  // 入射光线方向（单位向量）
    lie_scalar_t p_B[3]; // 镜面参考点（本体系）
    lie_scalar_t n_B[3]; // 镜面法向（本体系，单位向量）
} ray_mirror_params_t;

// 优化参数结构
typedef struct {
    ray_mirror_params_t* rays;      // 光路参数数组
    lie_scalar_t* h_meas;           // 测量光程值
    int num_rays;                   // 光路数量
} mirror_optimization_args_t;

// 计算单条光路的光程
lie_scalar_t cal_light_dis(const lie_scalar_t* R, const lie_scalar_t* t, 
                          const lie_scalar_t* a, const lie_scalar_t* u,
                          const lie_scalar_t* p_B, const lie_scalar_t* n_B);

// 计算单条光路的雅可比矩阵（对SE(3)李代数的导数）
void calc_jacobian(const lie_scalar_t* R, const lie_scalar_t* t,
                   const lie_scalar_t* a, const lie_scalar_t* u,
                   const lie_scalar_t* p_B, const lie_scalar_t* n_B,
                   lie_scalar_t* jac); // 输出6维雅可比向量

// 残差函数：用于LM优化
void Get_RS2_IFM_MS_ResFcn(const opt_scalar_t* xi, void* args, opt_scalar_t* residual);

// 雅可比函数：用于LM优化  
void Get_RS2_IFM_MS_Jacobian(const opt_scalar_t* xi, void* args, opt_scalar_t* jacobian);

// 主函数：执行镜面反射位姿优化
int mirror_reflection_pose_optimization(
    const ray_mirror_params_t* rays,
    const lie_scalar_t* h_meas,
    int num_rays,
    const lie_scalar_t* xi0,        // 初始猜测 [tx, ty, tz, phix, phiy, phiz]
    lie_scalar_t* xi_result,        // 输出结果
    lie_scalar_t* final_residual,   // 最终残差
    int max_iter,
    lie_scalar_t tol_x,
    lie_scalar_t tol_r,
    int* iters                      // 输出实际迭代次数
);

#ifdef __cplusplus
}
#endif

#endif
