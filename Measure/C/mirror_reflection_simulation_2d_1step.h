#ifndef __MIRROR_REFLECTION_SIMULATION_2D_1STEP_H__
#define __MIRROR_REFLECTION_SIMULATION_2D_1STEP_H__

#ifdef __cplusplus
extern "C" {
#endif

#define _WIN32

// 2D镜面反射仿真参数结构体
typedef struct {
    double phi[3];      // phi参数
    double psi[3];      // psi参数
    double a_x[3];      // a的x坐标
    double a_y[3];      // a的y坐标
    double b[3];        // b参数
    double lambda;
    
    // 中间计算变量
    double gamma[3];
    double zeta[3];
    double dx[3];
    double dy[3];

    double grad_theta[3];

    double h_pred[3];
    double h_meas[3];

    double UTU[4];
    double UTU_inv[6];
    double v[3];

    double x, y, theta;
} MirrorSim2D_Params;

// 函数声明

void module_init(MirrorSim2D_Params* in_param);
void precal(double *h_meas);
void cal_3dof(double *h_meas);

#ifdef __cplusplus
}
#endif

#endif // __MIRROR_REFLECTION_SIMULATION_2D_1STEP_H__
