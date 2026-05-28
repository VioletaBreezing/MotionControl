#ifndef __MIRROR_REFLECTION_SIMULATION_2D_H__
#define __MIRROR_REFLECTION_SIMULATION_2D_H__

#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

// 定义标量类型
#ifndef LIE_SOPHUS_USE_FLOAT
typedef double lie_scalar_t;
#else
typedef float lie_scalar_t;
#endif

// 数学函数适配
#ifdef LIE_SOPHUS_USE_FLOAT
#define lie_fabs fabsf
#define lie_sqrt sqrtf
#define lie_sin my_sin
#define lie_cos my_sin
#define lie_atan my_atan
#define lie_mean(a, b) (((a) + (b)) * 0.5f)
#else
#define lie_fabs fabs
#define lie_sqrt sqrt
#define lie_sin my_sin
#define lie_cos my_cos
#define lie_atan my_atan
#define lie_mean(a, b) (((a) + (b)) * 0.5)
#endif

// 2D镜面反射仿真参数结构体
typedef struct {
    lie_scalar_t phi[3];      // phi参数
    lie_scalar_t psi[3];      // psi参数
    lie_scalar_t a_x[3];      // a的x坐标
    lie_scalar_t a_y[3];      // a的y坐标
    lie_scalar_t b[3];        // b参数
    
    // 中间计算变量
    lie_scalar_t gamma[3];
    lie_scalar_t zeta[3];
    lie_scalar_t beta[3];
    lie_scalar_t sin_zeta[3];
    lie_scalar_t cos_zeta[3];
    lie_scalar_t f[3];
    lie_scalar_t dx[3];
    lie_scalar_t dy[3];
    lie_scalar_t g[3];
} MirrorSim2D_Params;

// 函数声明

/**
 * @brief 计算光程距离 h = f .* g
 * @param params 参数结构体指针
 * @param h 输出结果 [h1, h2, h3]
 */
// void cal_light_distance(const MirrorSim2D_Params* params, lie_scalar_t* h);

/**
 * @brief 更新参数结构体中的中间变量
 * @param x x方向位移
 * @param y y方向位移
 * @param rz 旋转角度
 * @param params 参数结构体指针（会被更新）
 */
// void update_params(lie_scalar_t x, lie_scalar_t y, lie_scalar_t rz, MirrorSim2D_Params* params);

// /**
//  * @brief 计算梯度
//  * @param params 参数结构体指针
//  * @param grad 输出梯度向量 [grad1, grad2, grad3]
//  */
// void cal_gradient(const MirrorSim2D_Params* params, lie_scalar_t* grad);

/**
 * @brief 计算3自由度位姿 (x, y, rz)
 * @param h_measure 测量的光程值 [h1, h2, h3]
 * @param params 参数结构体指针
 * @param rz0 初始rz猜测值
 * @param flag 标志位：0-需要估计rz，1-rz已知
 * @param tx 输出的x位移
 * @param ty 输出的y位移
 * @param rz_out 输出的旋转角度
 */
// void cal_3dof(const lie_scalar_t* h_measure, const MirrorSim2D_Params* params, 
//               lie_scalar_t rz0, int flag, 
//               lie_scalar_t* tx, lie_scalar_t* ty, lie_scalar_t* rz_out);

/**
 * @brief 使用梯度下降法计算3自由度位姿
 * @param h 测量的光程值 [h1, h2, h3]
 * @param params 参数结构体指针
 * @param rz0 初始rz猜测值
 * @param tx 输出的x位移
 * @param ty 输出的y位移
 * @param rz 输出的旋转角度
 * @param t 时间统计数组 [8个计时点]
 */
// void cal_3dof_grad(const lie_scalar_t* h, MirrorSim2D_Params* params, lie_scalar_t rz0,
//                    lie_scalar_t* tx, lie_scalar_t* ty, lie_scalar_t* rz, lie_scalar_t* t);

/**
 * @brief 运行主测试函数
 * @param N 迭代次数
 */
void mirror_sim_2d_main(int N);

#ifdef __cplusplus
}
#endif

#endif // __MIRROR_REFLECTION_SIMULATION_2D_H__
