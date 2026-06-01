#include "mirror_reflection_simulation_2d_1step.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// 添加max宏定义
#ifndef max
#define max(a,b) ((a) > (b) ? (a) : (b))
#endif

#ifdef _WIN32
#include <time.h>
#include <windows.h>
typedef LARGE_INTEGER timer_t;
timer_t freq;
#define timer_start(x)  do{QueryPerformanceCounter(&x);}while(0)
#define timer_stop(x)   do{QueryPerformanceCounter(&x);}while(0)
#define cal_time_us(x, y, t)  do{t = (double)(y.QuadPart - x.QuadPart) * 1000000.0 / freq.QuadPart;}\
                                while(0);
#define random_seed (unsigned int)time(NULL)
#elif defined(SOC_6678)
#include "BSP_int.h"
#include "BSP_timer.h"
typedef double timer_t;
#define timer_start
#define timer_stop
#define random_seed 0
#endif

#ifdef SOC_C6678
#pragma CODE_SECTION(cal_light_distance, ".sa_code")
#pragma CODE_SECTION(update_params, ".sa_code")
#pragma CODE_SECTION(cal_gradient, ".sa_code")
#pragma CODE_SECTION(cal_3dof, ".sa_code")
#pragma CODE_SECTION(cal_3dof_grad, ".sa_code")
#pragma CODE_SECTION(mirror_sim_2d_main, ".sa_code")
#endif

static const double PI = 3.14159265358979323846;

static MirrorSim2D_Params param;

void module_init(MirrorSim2D_Params* in_param)
{
    memcpy(&param, in_param, sizeof(MirrorSim2D_Params));
};

void cal_light_dis(double x, double y, double rz, 
    double *psi, double *phi, double *a_x, double *a_y, double *b,
    double *h)
{
    int i;
    for (i = 0; i < 3; i++)
    {
        double gamma = rz + psi[i] - phi[i];
        double beta = -cos(gamma);
        double f = beta / (2*beta*beta - 1);

        double zeta = (i == 0) ? (rz + psi[i] - PI/2) : (rz + psi[i]);
        double dx = x - a_x[i];
        double dy = y - a_y[i];
        double g = b[i] + dx * sin(zeta) - dy * cos(zeta);

        h[i] = f * g;
    }
}

// 梯度完全由上一时刻的值决定，可以放在预计算部分
void precal(double *h_meas)
{
    int i;
    for (i = 0; i < 3; i++) {
        param.dx[i] = param.x - param.a_x[i];
        param.dy[i] = param.y - param.a_y[i];
        param.h_meas[i] = h_meas[i];
    }

    // param.h_pred[0] = param.dx[0] + param.dy[0] * param.zeta[0] - param.b[0];
    param.grad_theta[0] = param.dx[0] * (param.zeta[0] - 9.0 * param.gamma[0]) + \
                          param.dy[0] * (1.0 + 9 * param.gamma[0] * param.zeta[0]) - \
                          9.0 * param.b[0] * param.gamma[0];
    for (i = 1; i < 3; i++) {
        // param.h_pred[i] = -param.dx[i] * param.zeta[i] + param.dy[i] - param.b[i];
        param.grad_theta[i] = -param.dx[i] * (1.0 + 9 * param.gamma[i] * param.zeta[i]) - \
                           param.dy[i] * (param.zeta[i] - 9.0 * param.gamma[i]) - \
                           9.0 * param.b[i] * param.gamma[i];
    }
    cal_light_dis(param.x, param.y, param.theta, param.psi, param.phi, param.a_x, param.a_y, param.b, param.h_pred);
    double a = 0;
}

void cal_3dof(double *h_meas, double *x, double *y, double *theta)
{
    // Update theta
    double residual[3];
    residual[0] = param.h_pred[0] - h_meas[0];
    residual[1] = param.h_pred[1] - h_meas[1];
    residual[2] = param.h_pred[2] - h_meas[2];

    double dtheta = residual[0] * param.grad_theta[0] + \
                    residual[1] * param.grad_theta[1] + \
                    residual[2] * param.grad_theta[2];
    param.theta -= param.lambda * dtheta;

    // Update angles
    param.zeta[0] = param.theta + param.psi[0];
    param.zeta[1] = param.theta + param.psi[1];
    param.zeta[2] = param.theta + param.psi[2];
    param.gamma[0] = param.zeta[0] - param.phi[0];
    param.gamma[1] = param.zeta[1] - param.phi[1];
    param.gamma[2] = param.zeta[2] - param.phi[2];

    // Update x, y
    param.UTU[0] = 1.0 + param.zeta[1] * param.zeta[1] + param.zeta[2] * param.zeta[2];
    param.UTU[1] = param.zeta[0] - param.zeta[1] - param.zeta[2];
    param.UTU[2] = param.UTU[1];
    param.UTU[3] = 2.0 + param.zeta[0] * param.zeta[0];

    double det = param.UTU[0] * param.UTU[3] - param.UTU[1] * param.UTU[2];
    param.UTU_inv[0] =  param.UTU[3] / det;
    param.UTU_inv[1] = -param.UTU[1] / det;
    param.UTU_inv[2] =  param.UTU_inv[1];
    param.UTU_inv[3] =  param.UTU[0] / det;

    param.v[0] = -param.a_x[0] - param.a_y[0] * param.zeta[0] -param.b[0];
    param.v[1] =  param.a_x[1] * param.zeta[1] - param.a_y[1] -param.b[1];
    param.v[2] =  param.a_x[2] * param.zeta[2] - param.a_y[2] -param.b[2];

    double c[3];
    c[0] = h_meas[0] - param.v[0];
    c[1] = h_meas[1] - param.v[1];
    c[2] = h_meas[2] - param.v[2];

    double d[2];
    d[0] = c[0] - param.zeta[1] * c[1] - param.zeta[2] * c[2];
    d[1] = param.zeta[0] * c[0] + c[1] + c[2];

    param.x = param.UTU_inv[0] * d[0] + param.UTU_inv[1] * d[1];
    param.y = param.UTU_inv[2] * d[0] + param.UTU_inv[3] * d[1];

    *x = param.x;
    *y = param.y;
    *theta = param.theta;
}

void self_test()
{
    srand(random_seed);

    MirrorSim2D_Params init_param;
    init_param.psi[0] = 2e-3 * 2.0 * (((double)rand() / RAND_MAX) - 0.5);
    init_param.psi[1] = 2e-3 * 2.0 * (((double)rand() / RAND_MAX) - 0.5);
    init_param.psi[2] = 2e-3 * 2.0 * (((double)rand() / RAND_MAX) - 0.5);

    init_param.phi[0] = 2e-3 * 2.0 * (((double)rand() / RAND_MAX) - 0.5);
    init_param.phi[1] = 2e-3 * 2.0 * (((double)rand() / RAND_MAX) - 0.5);
    init_param.phi[2] = 2e-3 * 2.0 * (((double)rand() / RAND_MAX) - 0.5);

    // a = {[-0.2; 0.1], [-0.1; -0.2], [0.1; -0.2]}
    init_param.a_x[0] = -0.2; init_param.a_y[0] =  0.1;
    init_param.a_x[1] = -0.1; init_param.a_y[1] = -0.2;
    init_param.a_x[2] =  0.1; init_param.a_y[2] = -0.2;
    
    // b = [0.1; 0.1; 0.1]
    init_param.b[0] = 0.1;
    init_param.b[1] = 0.1;
    init_param.b[2] = 0.1;

    init_param.lambda = 1;

    /* 生成期望轨迹 */
    double Ts = 2e-4;
    int trajectory_len = (int)(50/Ts);
    double *traj_x = (double*)malloc(trajectory_len * sizeof(double));
    double *traj_y = (double*)malloc(trajectory_len * sizeof(double));
    double *traj_theta = (double*)malloc(trajectory_len * sizeof(double));

    if (traj_x == NULL || traj_y == NULL || traj_theta == NULL) {
        printf("Memory allocation failed!\n");
        return;
    }

    int i;
    for(i = 0; i < trajectory_len; i++) {
        traj_x[i] = 1e-3 * cos(2 * PI * i * Ts);
        traj_y[i] = 1e-1 * sin(2 * PI * i * Ts);
        traj_theta[i] = 1e-3 * sin(2 * PI * i * Ts);
    }

    double h[3];
    cal_light_dis(traj_x[0], traj_y[0], traj_theta[0], 
                init_param.psi, init_param.phi, init_param.a_x, init_param.a_y,
                init_param.b, h);
    init_param.x = traj_x[0];
    init_param.y = traj_y[0];
    init_param.theta = traj_theta[0];

    module_init(&init_param);
    precal(h);

    /* 性能评估 */
#ifdef _WIN32
    QueryPerformanceFrequency(&freq);
#endif
    timer_t start_tic, end_est_tic, end_precal_tic;
    double max_est_time = -1, max_precal_time = -1, avg_time, used_est_time, used_precal_time, total_time = 0;
    double err[3]; 
    double max_err[3] = {-1, -1, -1};

    double x_est, y_est, theta_est;

    /* 开始仿真 */
    for (i = 0; i < trajectory_len; i++) {
        // 计算理论光程
        cal_light_dis(traj_x[i], traj_y[i], traj_theta[i], 
                init_param.psi, init_param.phi, init_param.a_x, init_param.a_y,
                init_param.b, h);

        timer_start(start_tic);
        cal_3dof(h, &x_est, &y_est, &theta_est);
        timer_stop(end_est_tic);
        cal_time_us(start_tic, end_est_tic, used_est_time);

        precal(h);
        timer_stop(end_precal_tic);
        cal_time_us(start_tic, end_precal_tic, used_precal_time);

        err[0] = fabs(x_est - traj_x[i]);
        err[1] = fabs(y_est - traj_y[i]);
        err[2] = fabs(theta_est - traj_theta[i]);

        max_err[0] = max(max_err[0], err[0]);
        max_err[1] = max(max_err[1], err[1]);
        max_err[2] = max(max_err[2], err[2]);

        max_est_time = max(max_est_time, used_est_time);
        max_precal_time = max(max_precal_time, used_precal_time);
        total_time += used_precal_time;
    }

    avg_time = total_time / trajectory_len;
    printf("max_est_time: %.3f us\n", max_est_time);
    printf("max_precal_time: %.3f us\n", max_precal_time);
    printf("avg_time: %.3f us\n", avg_time);
    printf("max_err: %.2e %.2e %.2e nm(nrad)\n", max_err[0]*1e9, max_err[1]*1e9, max_err[2]*1e9);

}



