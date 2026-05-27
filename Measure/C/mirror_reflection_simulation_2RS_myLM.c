#include "mirror_reflection_simulation_2RS_myLM.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// 数学函数适配
#ifdef LIE_SOPHUS_USE_FLOAT
#define lie_fabs fabsf
#define lie_sqrt sqrtf
#else
#define lie_fabs fabs
#define lie_sqrt sqrt
#endif

#ifdef SOC_C6678
#pragma CODE_SECTION(cal_light_dis, ".sa_code")
#pragma CODE_SECTION(skew_matrix, ".sa_code")
#pragma CODE_SECTION(calc_jacobian, ".sa_code")
#pragma CODE_SECTION(Get_RS2_IFM_MS_ResFcn, ".sa_code")
#pragma CODE_SECTION(Get_RS2_IFM_MS_Jacobian, ".sa_code")
#pragma CODE_SECTION(mirror_reflection_pose_optimization, ".sa_code")
#endif

// 计算单条光路的光程
lie_scalar_t cal_light_dis(const lie_scalar_t* R, const lie_scalar_t* t, 
                          const lie_scalar_t* a, const lie_scalar_t* u,
                          const lie_scalar_t* p_B, const lie_scalar_t* n_B) {
    // N = R * n_B
    lie_scalar_t N[3];
    N[0] = R[0]*n_B[0] + R[1]*n_B[1] + R[2]*n_B[2];
    N[1] = R[3]*n_B[0] + R[4]*n_B[1] + R[5]*n_B[2];
    N[2] = R[6]*n_B[0] + R[7]*n_B[1] + R[8]*n_B[2];
    
    // beta = dot(u, N)
    lie_scalar_t beta = u[0]*N[0] + u[1]*N[1] + u[2]*N[2];
    
    // ga = dot(p_B, n_B)
    lie_scalar_t ga = p_B[0]*n_B[0] + p_B[1]*n_B[1] + p_B[2]*n_B[2];
    
    // d = t - a
    lie_scalar_t d[3] = {t[0] - a[0], t[1] - a[1], t[2] - a[2]};
    
    // de = dot(d, N)
    lie_scalar_t de = d[0]*N[0] + d[1]*N[1] + d[2]*N[2];
    
    // f = beta / (2*beta^2 - 1)
    lie_scalar_t beta_sq = beta * beta;
    lie_scalar_t denom = 2.0 * beta_sq - 1.0;
    if (lie_fabs(denom) < 1e-12) {
        return 0.0; // 避免除零
    }
    lie_scalar_t f = beta / denom;
    
    // g = ga + de
    lie_scalar_t g = ga + de;
    
    // h = f * g
    return f * g;
}

// 计算反对称矩阵
void skew_matrix(const lie_scalar_t* a, lie_scalar_t* skew) {
    skew[0] = 0.0;      skew[1] = -a[2];   skew[2] = a[1];
    skew[3] = a[2];     skew[4] = 0.0;     skew[5] = -a[0];
    skew[6] = -a[1];    skew[7] = a[0];    skew[8] = 0.0;
}

// 计算单条光路的雅可比矩阵（对SE(3)李代数的导数）
void calc_jacobian(const lie_scalar_t* R, const lie_scalar_t* t,
                   const lie_scalar_t* a, const lie_scalar_t* u,
                   const lie_scalar_t* p_B, const lie_scalar_t* n_B,
                   lie_scalar_t* jac) { // 输出6维雅可比向量 [dt_x, dt_y, dt_z, dphi_x, dphi_y, dphi_z]
    
    // N = R * n_B
    lie_scalar_t N[3];
    N[0] = R[0]*n_B[0] + R[1]*n_B[1] + R[2]*n_B[2];
    N[1] = R[3]*n_B[0] + R[4]*n_B[1] + R[5]*n_B[2];
    N[2] = R[6]*n_B[0] + R[7]*n_B[1] + R[8]*n_B[2];
    
    // beta = dot(u, N)
    lie_scalar_t beta = u[0]*N[0] + u[1]*N[1] + u[2]*N[2];
    
    // ga = dot(p_B, n_B)
    lie_scalar_t ga = p_B[0]*n_B[0] + p_B[1]*n_B[1] + p_B[2]*n_B[2];
    
    // d = t - a
    lie_scalar_t d[3] = {t[0] - a[0], t[1] - a[1], t[2] - a[2]};
    
    // de = dot(d, N)
    lie_scalar_t de = d[0]*N[0] + d[1]*N[1] + d[2]*N[2];
    
    // f = beta / (2*beta^2 - 1)
    lie_scalar_t beta_sq = beta * beta;
    lie_scalar_t denom = 2.0 * beta_sq - 1.0;
    if (lie_fabs(denom) < 1e-12) {
        // 初始化为0
        int i;
        for (i = 0; i < 6; i++) {
            jac[i] = 0.0;
        }
        return;
    }
    lie_scalar_t f = beta / denom;
    
    // g = ga + de
    lie_scalar_t g = ga + de;
    
    // alph = -(2*beta^2 + 1) / (2*beta^2 - 1)^2
    lie_scalar_t alph = -(2.0 * beta_sq + 1.0) / (denom * denom);
    
    // pt = f * N' (平移部分的雅可比)
    jac[0] = f * N[0];  // dt_x
    jac[1] = f * N[1];  // dt_y  
    jac[2] = f * N[2];  // dt_z
    
    // 计算旋转部分的雅可比
    // ph = -transpose(alph * g * u + f * d) * skew_matrix(N)
    lie_scalar_t temp_vec[3];
    temp_vec[0] = alph * g * u[0] + f * d[0];
    temp_vec[1] = alph * g * u[1] + f * d[1];
    temp_vec[2] = alph * g * u[2] + f * d[2];
    
    // 计算 skew_matrix(N)
    lie_scalar_t skew_N[9];
    skew_matrix(N, skew_N);
    
    // ph = -temp_vec' * skew_N
    // 这等价于 ph = -(skew_N' * temp_vec)
    // 但更简单的方式是直接计算
    jac[3] = -(temp_vec[0] * skew_N[0] + temp_vec[1] * skew_N[3] + temp_vec[2] * skew_N[6]); // dphi_x
    jac[4] = -(temp_vec[0] * skew_N[1] + temp_vec[1] * skew_N[4] + temp_vec[2] * skew_N[7]); // dphi_y
    jac[5] = -(temp_vec[0] * skew_N[2] + temp_vec[1] * skew_N[5] + temp_vec[2] * skew_N[8]); // dphi_z
}

// 残差函数：用于LM优化
void Get_RS2_IFM_MS_ResFcn(const opt_scalar_t* xi, void* args, opt_scalar_t* residual) {
    mirror_optimization_args_t* opt_args = (mirror_optimization_args_t*)args;
    int num_rays = opt_args->num_rays;
    
    // 将xi转换为SE(3)变换矩阵
    lie_scalar_t T[16];
    se3_to_SE3((const lie_scalar_t*)xi, T);
    
    lie_scalar_t R[9] = {T[0], T[1], T[2], T[4], T[5], T[6], T[8], T[9], T[10]};
    lie_scalar_t t[3] = {T[3], T[7], T[11]};
    
    // 计算每条光路的预测光程
    // int i;
    // for (i = 0; i < num_rays; i++) {
    //     lie_scalar_t h_pred = cal_light_dis(R, t, 
    //                                       opt_args->rays[i].a, opt_args->rays[i].u,
    //                                       opt_args->rays[i].p_B, opt_args->rays[i].n_B);
    //     residual[i] = h_pred - opt_args->h_meas[i];
    // }
    lie_scalar_t h_pred = cal_light_dis(R, t, opt_args->rays[0].a, opt_args->rays[0].u, opt_args->rays[0].p_B, opt_args->rays[0].n_B);
    residual[0] = h_pred - opt_args->h_meas[0];
    h_pred = cal_light_dis(R, t, opt_args->rays[1].a, opt_args->rays[1].u, opt_args->rays[1].p_B, opt_args->rays[1].n_B);
    residual[1] = h_pred - opt_args->h_meas[1];
    h_pred = cal_light_dis(R, t, opt_args->rays[2].a, opt_args->rays[2].u, opt_args->rays[2].p_B, opt_args->rays[2].n_B);
    residual[2] = h_pred - opt_args->h_meas[2];
}

// 雅可比函数：用于LM优化
void Get_RS2_IFM_MS_Jacobian(const opt_scalar_t* xi, void* args, opt_scalar_t* jacobian) {
    mirror_optimization_args_t* opt_args = (mirror_optimization_args_t*)args;
    int num_rays = opt_args->num_rays;
    
    // 将xi转换为SE(3)变换矩阵
    lie_scalar_t T[16];
    se3_to_SE3((const lie_scalar_t*)xi, T);
    
    lie_scalar_t R[9] = {T[0], T[1], T[2], T[4], T[5], T[6], T[8], T[9], T[10]};
    lie_scalar_t t[3] = {T[3], T[7], T[11]};
    
    // 计算每条光路的雅可比
    int i;
    for (i = 0; i < num_rays; i++) {
        lie_scalar_t jac_6d[6];
        calc_jacobian(R, t,
                     opt_args->rays[i].a, opt_args->rays[i].u,
                     opt_args->rays[i].p_B, opt_args->rays[i].n_B,
                     jac_6d);
        
        // 将6D雅可比复制到输出矩阵中
        jacobian[i * 6 + 0] = jac_6d[0];
        jacobian[i * 6 + 1] = jac_6d[1];
        jacobian[i * 6 + 2] = 0.;
        jacobian[i * 6 + 3] = 0.;
        jacobian[i * 6 + 4] = 0.;
        jacobian[i * 6 + 5] = jac_6d[5];
    }
}

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
    int *iters) {
    
    // 准备优化参数
    mirror_optimization_args_t opt_args;
    opt_args.rays = (ray_mirror_params_t*)rays;
    opt_args.h_meas = (lie_scalar_t*)h_meas;
    opt_args.num_rays = num_rays;
    
    // 设置LM参数
    opt_scalar_t x_result[6];
    opt_scalar_t r_result[6];
    lm_result_t result;
    result.x_hat = x_result;
    result.resnorm = r_result;
    
    int status = Optimizer_LM(
        Get_RS2_IFM_MS_ResFcn,
        &opt_args,
        Get_RS2_IFM_MS_Jacobian,
        &opt_args,
        (const opt_scalar_t*)xi0,
        6,          // param_size: 6DoF位姿
        num_rays,   // residual_size: 光路数量
        tol_x,
        tol_r,
        1e10,       // TolLamda
        max_iter,
        0,          // debug: 关闭
        &result
    );

    *iters = result.iter;
    
    if (status == 0) {
        // 成功优化
        memcpy(xi_result, result.x_hat, 6 * sizeof(lie_scalar_t));
        if (final_residual != NULL) {
            memcpy(final_residual, result.resnorm, num_rays * sizeof(lie_scalar_t));
        }
    }
    
    // 清理内存（动态分配模式下）
#if !defined(OPT_LM_FORBIDDEN_DYNAMIC_ALLOCATION) || OPT_LM_FORBIDDEN_DYNAMIC_ALLOCATION == 0
    free(result.x_hat);
    free(result.resnorm);
#endif
    
    return status;
}
