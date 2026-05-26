#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#ifndef min
#define min(a,b) (((a)<(b))?(a):(b))
#endif  

#ifndef max
#define max(a,b) (((a)>(b))?(a):(b))
#endif

#ifdef _SIMULAE_IN_X86_

#include "Measure/C/mirror_reflection_simulation_2RS_myLM.h"
#include <windows.h>
#include <time.h>

int main() {
    printf("=== Testing Mirror Reflection Pose Optimization with Real Initial Guess ===\n");

    srand((unsigned int)time(NULL));
    
    // 设置真实位姿参数（与MATLAB一致）
    lie_scalar_t theta_x = 0.0;           // 绕X轴旋转（设为0）
    lie_scalar_t theta_y = 0.0;           // 绕Y轴旋转（设为0）
    lie_scalar_t theta_z = 20e-6;         // 绕Z轴旋转（非零）
    
    lie_scalar_t tx = 1e-3;               // X平移
    lie_scalar_t ty = -1e-3;              // Y平移  
    lie_scalar_t tz = 0.0;                // Z平移（设为0）
    
    // 构建真实旋转矩阵 Rz * Ry * Rx
    lie_scalar_t cos_x = cos(theta_x), sin_x = sin(theta_x);
    lie_scalar_t cos_y = cos(theta_y), sin_y = sin(theta_y);
    lie_scalar_t cos_z = cos(theta_z), sin_z = sin(theta_z);
    
    lie_scalar_t Rx[9] = {1, 0, 0, 0, cos_x, -sin_x, 0, sin_x, cos_x};
    lie_scalar_t Ry[9] = {cos_y, 0, sin_y, 0, 1, 0, -sin_y, 0, cos_y};
    lie_scalar_t Rz[9] = {cos_z, -sin_z, 0, sin_z, cos_z, 0, 0, 0, 1};
    
    // R = Rz * Ry * Rx
    lie_scalar_t R_temp[9], R_true[9];
    MatrixMath_Multiply(Rz, Ry, 3, 3, 3, R_temp);
    MatrixMath_Multiply(R_temp, Rx, 3, 3, 3, R_true);
    
    lie_scalar_t t_true[3] = {tx, ty, tz};
    
    // 刚体与镜面定义
    lie_scalar_t Lx = 0.2, Ly = 0.2, Lz = 0.1;
    lie_scalar_t d = 0.1;
    
    // 镜面参考点（本体系）
    lie_scalar_t p0_B[3] = {-Lx/2, 0, 0};
    lie_scalar_t p1_B[3] = {-d/2, -Ly/2, 0};
    lie_scalar_t p2_B[3] = {d/2, -Ly/2, 0};
    
    // 镜面法向（本体系）
    lie_scalar_t n0_B[3] = {-1, 0.1, -0.1};
    lie_scalar_t n1_B[3] = {0.05, -1, -0.02};
    lie_scalar_t n2_B[3] = {-0.12, -1, 0.08};
    
    // 归一化法向
    lie_scalar_t norm_n0 = sqrt(n0_B[0]*n0_B[0] + n0_B[1]*n0_B[1] + n0_B[2]*n0_B[2]);
    lie_scalar_t norm_n1 = sqrt(n1_B[0]*n1_B[0] + n1_B[1]*n1_B[1] + n1_B[2]*n1_B[2]);
    lie_scalar_t norm_n2 = sqrt(n2_B[0]*n2_B[0] + n2_B[1]*n2_B[1] + n2_B[2]*n2_B[2]);
    n0_B[0] /= norm_n0; n0_B[1] /= norm_n0; n0_B[2] /= norm_n0;
    n1_B[0] /= norm_n1; n1_B[1] /= norm_n1; n1_B[2] /= norm_n1;
    n2_B[0] /= norm_n2; n2_B[1] /= norm_n2; n2_B[2] /= norm_n2;
    
    // 入射光线
    lie_scalar_t a0[3] = {-0.3, 0.02, 0.0};
    lie_scalar_t u0[3] = {1.0, 0.01, -0.01};
    lie_scalar_t a1[3] = {-d/2, -0.3, 0.0};
    lie_scalar_t u1[3] = {-0.01, 1.0, 0.02};
    lie_scalar_t a2[3] = {d/2, -0.3, 0.0};
    lie_scalar_t u2[3] = {0.04, 1.0, 0.03};
    
    // 归一化入射方向
    lie_scalar_t norm_u0 = sqrt(u0[0]*u0[0] + u0[1]*u0[1] + u0[2]*u0[2]);
    lie_scalar_t norm_u1 = sqrt(u1[0]*u1[0] + u1[1]*u1[1] + u1[2]*u1[2]);
    lie_scalar_t norm_u2 = sqrt(u2[0]*u2[0] + u2[1]*u2[1] + u2[2]*u2[2]);
    u0[0] /= norm_u0; u0[1] /= norm_u0; u0[2] /= norm_u0;
    u1[0] /= norm_u1; u1[1] /= norm_u1; u1[2] /= norm_u1;
    u2[0] /= norm_u2; u2[1] /= norm_u2; u2[2] /= norm_u2;
    
    // 定义光路参数
    ray_mirror_params_t rays[3];
    
    // 光路0
    memcpy(rays[0].a, a0, 3 * sizeof(lie_scalar_t));
    memcpy(rays[0].u, u0, 3 * sizeof(lie_scalar_t));
    memcpy(rays[0].p_B, p0_B, 3 * sizeof(lie_scalar_t));
    memcpy(rays[0].n_B, n0_B, 3 * sizeof(lie_scalar_t));
    
    // 光路1
    memcpy(rays[1].a, a1, 3 * sizeof(lie_scalar_t));
    memcpy(rays[1].u, u1, 3 * sizeof(lie_scalar_t));
    memcpy(rays[1].p_B, p1_B, 3 * sizeof(lie_scalar_t));
    memcpy(rays[1].n_B, n1_B, 3 * sizeof(lie_scalar_t));
    
    // 光路2
    memcpy(rays[2].a, a2, 3 * sizeof(lie_scalar_t));
    memcpy(rays[2].u, u2, 3 * sizeof(lie_scalar_t));
    memcpy(rays[2].p_B, p2_B, 3 * sizeof(lie_scalar_t));
    memcpy(rays[2].n_B, n2_B, 3 * sizeof(lie_scalar_t));
    
    // 计算真实光程值（使用真实位姿）
    lie_scalar_t h_meas[3];
    h_meas[0] = cal_light_dis(R_true, t_true, a0, u0, p0_B, n0_B);
    h_meas[1] = cal_light_dis(R_true, t_true, a1, u1, p1_B, n1_B);
    h_meas[2] = cal_light_dis(R_true, t_true, a2, u2, p2_B, n2_B);
    
    printf("True measurements:\n");
    printf("h0 = %.10f m\n", h_meas[0]);
    printf("h1 = %.10f m\n", h_meas[1]);
    printf("h2 = %.10f m\n", h_meas[2]);
    
    // 使用MATLAB中的初始猜测（故意设错，但不是零）
    lie_scalar_t xi0[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // [tx, ty, tz, phix, phiy, phiz]
    
    printf("\nInitial guess: [%.6f, %.6f, %.6f, %.6f, %.6f, %.6f]\n",
           xi0[0], xi0[1], xi0[2], xi0[3], xi0[4], xi0[5]);
    
    // 执行优化并计时 (使用高精度计时器)
    lie_scalar_t xi_result[6];
    lie_scalar_t final_residual[3];
    
    LARGE_INTEGER frequency, start_time, start_time_iter, end_time;
    double max_exec_time = 1e12;
    QueryPerformanceFrequency(&frequency);
    QueryPerformanceCounter(&start_time);

    int times_cnt = 20000;
    int status = 0;
    lie_scalar_t T[16];
    lie_scalar_t R_est[9];
    lie_scalar_t t_est[3];
    lie_scalar_t max_trans_error[3] = {0.0, 0.0, 0.0};
    lie_scalar_t max_rot_error[3] = {0.0, 0.0, 0.0};
    int i;
    for(i = 0; i < times_cnt && status == 0; i++)
    {
        QueryPerformanceCounter(&start_time_iter);

        theta_z = (lie_scalar_t)rand() / RAND_MAX * 2.0 - 1.0;
        theta_z *= 2e-3; // 旋转范围在±1 mrad内
        
        tx = (lie_scalar_t)rand() / RAND_MAX * 2.0 - 1.0;
        ty = (lie_scalar_t)rand() / RAND_MAX * 2.0 - 1.0;
        tx *= 13e-2; // 平移范围在±1 mm内
        ty *= 13e-2;

        h_meas[0] = cal_light_dis(R_true, t_true, a0, u0, p0_B, n0_B);
        h_meas[1] = cal_light_dis(R_true, t_true, a1, u1, p1_B, n1_B);
        h_meas[2] = cal_light_dis(R_true, t_true, a2, u2, p2_B, n2_B);

        cos_x = cos(theta_x), sin_x = sin(theta_x);
        cos_y = cos(theta_y), sin_y = sin(theta_y);
        cos_z = cos(theta_z), sin_z = sin(theta_z);
        
        Rx[0] = 1; Rx[1] = 0;     Rx[2] = 0;
        Rx[3] = 0; Rx[4] = cos_x; Rx[5] = -sin_x;
        Rx[6] = 0; Rx[7] = sin_x; Rx[8] = cos_x;

        Ry[0] = cos_y;  Ry[1] = 0; Ry[2] = sin_y;
        Ry[3] = 0;      Ry[4] = 1; Ry[5] = 0;
        Ry[6] = -sin_y; Ry[7] = 0; Ry[8] = cos_y;

        Rz[0] = cos_z; Rz[1] = -sin_z; Rz[2] = 0;
        Rz[3] = sin_z; Rz[4] = cos_z;  Rz[5] = 0;
        Rz[6] = 0;     Rz[7] = 0;      Rz[8] = 1;
        
        MatrixMath_Multiply(Rz, Ry, 3, 3, 3, R_temp);
        MatrixMath_Multiply(R_temp, Rx, 3, 3, 3, R_true);
    
        t_true[0] = tx; t_true[1] = ty; t_true[2] = tz;

        h_meas[0] = cal_light_dis(R_true, t_true, a0, u0, p0_B, n0_B);
        h_meas[1] = cal_light_dis(R_true, t_true, a1, u1, p1_B, n1_B);
        h_meas[2] = cal_light_dis(R_true, t_true, a2, u2, p2_B, n2_B);

        status = mirror_reflection_pose_optimization(
            rays,
            h_meas,
            3,
            xi0,
            xi_result,
            final_residual,
            100,     // max_iter
            1e-11,  // tol_x (更严格的收敛条件)
            1e-11   // tol_r (更严格的收敛条件)
        );

        se3_to_SE3(xi_result, T);
        
        R_est[0] = T[0]; R_est[1] = T[1]; R_est[2] = T[2];
        R_est[3] = T[4]; R_est[4] = T[5]; R_est[5] = T[6];
        R_est[6] = T[8]; R_est[7] = T[9]; R_est[8] = T[10];

        t_est[0] = T[3]; t_est[1] = T[7]; t_est[2] = T[11];

        QueryPerformanceCounter(&end_time);
        double iter_time_used = (double)(end_time.QuadPart - start_time_iter.QuadPart) / frequency.QuadPart;
        max_exec_time = min(max_exec_time, iter_time_used);

        max_trans_error[0] = max(max_trans_error[0], fabs(t_est[0] - tx));
        max_trans_error[1] = max(max_trans_error[1], fabs(t_est[1] - ty));
        max_trans_error[2] = max(max_trans_error[2], fabs(t_est[2] - tz));

        max_rot_error[0] = max(max_rot_error[0], fabs(xi_result[3] - theta_x));
        max_rot_error[1] = max(max_rot_error[1], fabs(xi_result[4] - theta_y));
        max_rot_error[2] = max(max_rot_error[2], fabs(xi_result[5] - theta_z));
    }
    
    QueryPerformanceCounter(&end_time);
    double cpu_time_used = (double)(end_time.QuadPart - start_time.QuadPart) / frequency.QuadPart / times_cnt;
    
    if (status == 0) {
        // lie_scalar_t T[16];
        // se3_to_SE3(xi_result, T);
        // lie_scalar_t R_est[9] = {T[0], T[1], T[2], T[4], T[5], T[6], T[8], T[9], T[10]};
        // lie_scalar_t t_est[3] = {T[3], T[7], T[11]};

        printf("\n✅ Optimization successful!\n");
        printf("Estimated pose:\n");
        printf("Translation (mm): [%.7f, %.7f, %.7f]\n", t_est[0] * 1e3, t_est[1] * 1e3, t_est[2] * 1e3);
        printf("Rotation (mrad):  [%.7f, %.7f, %.7f]\n", xi_result[3] * 1e3, xi_result[4] * 1e3, xi_result[5] * 1e3);
        
        printf("\nTrue pose:\n");
        printf("Translation (mm): [%.7f, %.7f, %.7f]\n", tx * 1e3, ty * 1e3, tz * 1e3);
        printf("Rotation (mrad):  [%.7f, %.7f, %.7f]\n", theta_x * 1e3, theta_y * 1e3, theta_z * 1e3);
        
        // 计算误差
        lie_scalar_t trans_error[3] = {t_est[0] - tx, t_est[1] - ty, t_est[2] - tz};
        lie_scalar_t rot_error[3] = {xi_result[3] - theta_x, xi_result[4] - theta_y, xi_result[5] - theta_z};
        
        printf("\nErrors:\n");
        printf("Translation error (nm): [%.2e, %.2e, %.2e]\n", 
               trans_error[0]*1e9, trans_error[1]*1e9, trans_error[2]*1e9);
        printf("Rotation error (nrad):  [%.2e, %.2e, %.2e]\n", 
               rot_error[0]*1e9, rot_error[1]*1e9, rot_error[2]*1e9);

        printf("\nMaximum translation error (nm): [%.3e, %.3e, %.3e]\n", 
               max_trans_error[0]*1e9, max_trans_error[1]*1e9, max_trans_error[2]*1e9);
        printf("Maximum rotation error (nrad):  [%.3e, %.3e, %.3e]\n", 
               max_rot_error[0]*1e9, max_rot_error[1]*1e9, max_rot_error[2]*1e9);
        
        printf("\nFinal residuals (nm):   [%.2e, %.2e, %.2e]\n", 
                final_residual[0]*1e9, final_residual[1]*1e9, final_residual[2]*1e9);
        
        // 检查是否成功收敛到真实值
        lie_scalar_t trans_error_norm = sqrt(trans_error[0]*trans_error[0] + trans_error[1]*trans_error[1] + trans_error[2]*trans_error[2]);
        lie_scalar_t rot_error_norm = sqrt(rot_error[0]*rot_error[0] + rot_error[1]*rot_error[1] + rot_error[2]*rot_error[2]);
        
        if (trans_error_norm < 1e-9 && rot_error_norm < 1e-8) {
            printf("\n🎉 Perfect convergence achieved!\n");
        } else if (trans_error_norm < 1e-6 && rot_error_norm < 1e-6) {
            printf("\n👍 Good convergence achieved!\n");
        } else {
            printf("\n⚠️  Convergence may be suboptimal. Consider adjusting parameters.\n");
        }
        
        printf("\n⏱️  Execution time: %.6f us\n", cpu_time_used * 1e6);
        printf("Maximum execution time: %.6f us\n", max_exec_time * 1e6);
    } else {
        printf("❌ Optimization failed with status: %d\n", status);
        printf("Execution time: %.6f seconds\n", cpu_time_used);
        printf("Status codes:\n");
        printf("  0 = Success\n");
        printf(" -1 = Invalid input\n");
        printf(" -2 = Memory allocation failed\n");
        printf(" -3 = Matrix inversion failed\n");
        printf(" -4 = Lambda too large\n");
        printf("\n⏱️  Execution time: %.6f us\n", cpu_time_used * 1e6);
    }
    
    return 0;
}

#else
#include "BSP_uart.h"
#include "BSP_eth.h"
#include "BSP_shm.h"
#include "BSP_int.h"
#include "BSP_timer.h"

#include "DRDBLG.h"
#include "DRDBMT.h"

#include "DRIPCM.h"
#include "DRIPCD.h"
#include "mirror_reflection_simulation_2RS_myLM.h"

#pragma CODE_SECTION(test_ifm, ".sa_code")

int test_ifm() {
    LOG_CONSOLE("=== Testing Mirror Reflection Pose Optimization with Real Initial Guess ===\n");

    srand(0);
    
    // 设置真实位姿参数（与MATLAB一致）
    lie_scalar_t theta_x = 0.0;           // 绕X轴旋转（设为0）
    lie_scalar_t theta_y = 0.0;           // 绕Y轴旋转（设为0）
    lie_scalar_t theta_z = 20e-6;         // 绕Z轴旋转（非零）
    
    lie_scalar_t tx = 1e-3;               // X平移
    lie_scalar_t ty = -1e-3;              // Y平移  
    lie_scalar_t tz = 0.0;                // Z平移（设为0）
    
    // 构建真实旋转矩阵 Rz * Ry * Rx
    lie_scalar_t cos_x = cos(theta_x), sin_x = sin(theta_x);
    lie_scalar_t cos_y = cos(theta_y), sin_y = sin(theta_y);
    lie_scalar_t cos_z = cos(theta_z), sin_z = sin(theta_z);
    
    lie_scalar_t Rx[9] = {1, 0, 0, 0, cos_x, -sin_x, 0, sin_x, cos_x};
    lie_scalar_t Ry[9] = {cos_y, 0, sin_y, 0, 1, 0, -sin_y, 0, cos_y};
    lie_scalar_t Rz[9] = {cos_z, -sin_z, 0, sin_z, cos_z, 0, 0, 0, 1};
    
    // R = Rz * Ry * Rx
    lie_scalar_t R_temp[9], R_true[9];
    MatrixMath_Multiply(Rz, Ry, 3, 3, 3, R_temp);
    MatrixMath_Multiply(R_temp, Rx, 3, 3, 3, R_true);
    
    lie_scalar_t t_true[3] = {tx, ty, tz};
    
    // 刚体与镜面定义
    lie_scalar_t Lx = 0.2, Ly = 0.2, Lz = 0.1;
    lie_scalar_t d = 0.1;
    
    // 镜面参考点（本体系）
    lie_scalar_t p0_B[3] = {-Lx/2, 0, 0};
    lie_scalar_t p1_B[3] = {-d/2, -Ly/2, 0};
    lie_scalar_t p2_B[3] = {d/2, -Ly/2, 0};
    
    // 镜面法向（本体系）
    lie_scalar_t n0_B[3] = {-1, 0.1, -0.1};
    lie_scalar_t n1_B[3] = {0.05, -1, -0.02};
    lie_scalar_t n2_B[3] = {-0.12, -1, 0.08};
    
    // 归一化法向
    lie_scalar_t norm_n0 = sqrt(n0_B[0]*n0_B[0] + n0_B[1]*n0_B[1] + n0_B[2]*n0_B[2]);
    lie_scalar_t norm_n1 = sqrt(n1_B[0]*n1_B[0] + n1_B[1]*n1_B[1] + n1_B[2]*n1_B[2]);
    lie_scalar_t norm_n2 = sqrt(n2_B[0]*n2_B[0] + n2_B[1]*n2_B[1] + n2_B[2]*n2_B[2]);
    n0_B[0] /= norm_n0; n0_B[1] /= norm_n0; n0_B[2] /= norm_n0;
    n1_B[0] /= norm_n1; n1_B[1] /= norm_n1; n1_B[2] /= norm_n1;
    n2_B[0] /= norm_n2; n2_B[1] /= norm_n2; n2_B[2] /= norm_n2;
    
    // 入射光线
    lie_scalar_t a0[3] = {-0.3, 0.02, 0.0};
    lie_scalar_t u0[3] = {1.0, 0.01, -0.01};
    lie_scalar_t a1[3] = {-d/2, -0.3, 0.0};
    lie_scalar_t u1[3] = {-0.01, 1.0, 0.02};
    lie_scalar_t a2[3] = {d/2, -0.3, 0.0};
    lie_scalar_t u2[3] = {0.04, 1.0, 0.03};
    
    // 归一化入射方向
    lie_scalar_t norm_u0 = sqrt(u0[0]*u0[0] + u0[1]*u0[1] + u0[2]*u0[2]);
    lie_scalar_t norm_u1 = sqrt(u1[0]*u1[0] + u1[1]*u1[1] + u1[2]*u1[2]);
    lie_scalar_t norm_u2 = sqrt(u2[0]*u2[0] + u2[1]*u2[1] + u2[2]*u2[2]);
    u0[0] /= norm_u0; u0[1] /= norm_u0; u0[2] /= norm_u0;
    u1[0] /= norm_u1; u1[1] /= norm_u1; u1[2] /= norm_u1;
    u2[0] /= norm_u2; u2[1] /= norm_u2; u2[2] /= norm_u2;
    
    // 定义光路参数
    ray_mirror_params_t rays[3];
    
    // 光路0
    memcpy(rays[0].a, a0, 3 * sizeof(lie_scalar_t));
    memcpy(rays[0].u, u0, 3 * sizeof(lie_scalar_t));
    memcpy(rays[0].p_B, p0_B, 3 * sizeof(lie_scalar_t));
    memcpy(rays[0].n_B, n0_B, 3 * sizeof(lie_scalar_t));
    
    // 光路1
    memcpy(rays[1].a, a1, 3 * sizeof(lie_scalar_t));
    memcpy(rays[1].u, u1, 3 * sizeof(lie_scalar_t));
    memcpy(rays[1].p_B, p1_B, 3 * sizeof(lie_scalar_t));
    memcpy(rays[1].n_B, n1_B, 3 * sizeof(lie_scalar_t));
    
    // 光路2
    memcpy(rays[2].a, a2, 3 * sizeof(lie_scalar_t));
    memcpy(rays[2].u, u2, 3 * sizeof(lie_scalar_t));
    memcpy(rays[2].p_B, p2_B, 3 * sizeof(lie_scalar_t));
    memcpy(rays[2].n_B, n2_B, 3 * sizeof(lie_scalar_t));
    
    // 计算真实光程值（使用真实位姿）
    lie_scalar_t h_meas[3];
    h_meas[0] = cal_light_dis(R_true, t_true, a0, u0, p0_B, n0_B);
    h_meas[1] = cal_light_dis(R_true, t_true, a1, u1, p1_B, n1_B);
    h_meas[2] = cal_light_dis(R_true, t_true, a2, u2, p2_B, n2_B);
    
    LOG_CONSOLE("True measurements:\n");
    LOG_CONSOLE("h0 = %.10f m\n", h_meas[0]);
    LOG_CONSOLE("h1 = %.10f m\n", h_meas[1]);
    LOG_CONSOLE("h2 = %.10f m\n", h_meas[2]);
    
    // 使用MATLAB中的初始猜测（故意设错，但不是零）
    lie_scalar_t xi0[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // [tx, ty, tz, phix, phiy, phiz]
    
    LOG_CONSOLE("\nInitial guess: [%.6f, %.6f, %.6f, %.6f, %.6f, %.6f]\n",
           xi0[0], xi0[1], xi0[2], xi0[3], xi0[4], xi0[5]);
    
    // 执行优化并计时 (使用高精度计时器)
    lie_scalar_t xi_result[6];
    lie_scalar_t final_residual[3];
    
    // LARGE_INTEGER frequency, start_time, start_time_iter, end_time;
    // double max_exec_time = 1e12;
    // QueryPerformanceFrequency(&frequency);
    // QueryPerformanceCounter(&start_time);
    int start_time, start_time_iter, end_time;
    double max_exec_time = 1e12;
    start_time = BSP_get_time();

    int times_cnt = 200;
    int status = 0;
    lie_scalar_t T[16];
    lie_scalar_t R_est[9];
    lie_scalar_t t_est[3];
    lie_scalar_t max_trans_error[3] = {0.0, 0.0, 0.0};
    lie_scalar_t max_rot_error[3] = {0.0, 0.0, 0.0};
    int i;
    for(i = 0; i < times_cnt && status == 0; i++)
    {
        start_time_iter = BSP_get_time();

        theta_z = (lie_scalar_t)rand() / RAND_MAX * 2.0 - 1.0;
        theta_z *= 2e-3; // 旋转范围在±1 mrad内
        
        tx = (lie_scalar_t)rand() / RAND_MAX * 2.0 - 1.0;
        ty = (lie_scalar_t)rand() / RAND_MAX * 2.0 - 1.0;
        tx *= 13e-2; // 平移范围在±1 mm内
        ty *= 13e-2;

        h_meas[0] = cal_light_dis(R_true, t_true, a0, u0, p0_B, n0_B);
        h_meas[1] = cal_light_dis(R_true, t_true, a1, u1, p1_B, n1_B);
        h_meas[2] = cal_light_dis(R_true, t_true, a2, u2, p2_B, n2_B);

        cos_x = cos(theta_x), sin_x = sin(theta_x);
        cos_y = cos(theta_y), sin_y = sin(theta_y);
        cos_z = cos(theta_z), sin_z = sin(theta_z);
        cos_x = max(-1, min(1, cos_x));
        sin_x = max(-1, min(1, sin_x));
        cos_y = max(-1, min(1, cos_y));
        sin_y = max(-1, min(1, sin_y));
        cos_z = max(-1, min(1, cos_z));
        sin_z = max(-1, min(1, sin_z));
        
        Rx[0] = 1; Rx[1] = 0;     Rx[2] = 0;
        Rx[3] = 0; Rx[4] = cos_x; Rx[5] = -sin_x;
        Rx[6] = 0; Rx[7] = sin_x; Rx[8] = cos_x;

        Ry[0] = cos_y;  Ry[1] = 0; Ry[2] = sin_y;
        Ry[3] = 0;      Ry[4] = 1; Ry[5] = 0;
        Ry[6] = -sin_y; Ry[7] = 0; Ry[8] = cos_y;

        Rz[0] = cos_z; Rz[1] = -sin_z; Rz[2] = 0;
        Rz[3] = sin_z; Rz[4] = cos_z;  Rz[5] = 0;
        Rz[6] = 0;     Rz[7] = 0;      Rz[8] = 1;
        
        MatrixMath_Multiply(Rz, Ry, 3, 3, 3, R_temp);
        MatrixMath_Multiply(R_temp, Rx, 3, 3, 3, R_true);
    
        t_true[0] = tx; t_true[1] = ty; t_true[2] = tz;

        h_meas[0] = cal_light_dis(R_true, t_true, a0, u0, p0_B, n0_B);
        h_meas[1] = cal_light_dis(R_true, t_true, a1, u1, p1_B, n1_B);
        h_meas[2] = cal_light_dis(R_true, t_true, a2, u2, p2_B, n2_B);

        status = mirror_reflection_pose_optimization(
            rays,
            h_meas,
            3,
            xi0,
            xi_result,
            final_residual,
            100,     // max_iter
            3e-10,  // tol_x (更严格的收敛条件)
            3e-10   // tol_r (更严格的收敛条件)
        );

        se3_to_SE3(xi_result, T);
        
        R_est[0] = T[0]; R_est[1] = T[1]; R_est[2] = T[2];
        R_est[3] = T[4]; R_est[4] = T[5]; R_est[5] = T[6];
        R_est[6] = T[8]; R_est[7] = T[9]; R_est[8] = T[10];

        t_est[0] = T[3]; t_est[1] = T[7]; t_est[2] = T[11];

        // QueryPerformanceCounter(&end_time);
        // double iter_time_used = (double)(end_time.QuadPart - start_time_iter.QuadPart) / frequency.QuadPart;
        // max_exec_time = min(max_exec_time, iter_time_used);
        end_time = BSP_get_time();
        double iter_time_used = BSP_calc_us_interval(start_time_iter, end_time);
        max_exec_time = min(max_exec_time, iter_time_used);

        max_trans_error[0] = max(max_trans_error[0], fabs(t_est[0] - tx));
        max_trans_error[1] = max(max_trans_error[1], fabs(t_est[1] - ty));
        max_trans_error[2] = max(max_trans_error[2], fabs(t_est[2] - tz));

        max_rot_error[0] = max(max_rot_error[0], fabs(xi_result[3] - theta_x));
        max_rot_error[1] = max(max_rot_error[1], fabs(xi_result[4] - theta_y));
        max_rot_error[2] = max(max_rot_error[2], fabs(xi_result[5] - theta_z));
    }
    
    // QueryPerformanceCounter(&end_time);
    // double cpu_time_used = (double)(end_time.QuadPart - start_time.QuadPart) / frequency.QuadPart / times_cnt;
    double cpu_time_used = BSP_calc_us_interval(start_time, end_time) / times_cnt;
    
    if (status == 0) {
        // lie_scalar_t T[16];
        // se3_to_SE3(xi_result, T);
        // lie_scalar_t R_est[9] = {T[0], T[1], T[2], T[4], T[5], T[6], T[8], T[9], T[10]};
        // lie_scalar_t t_est[3] = {T[3], T[7], T[11]};

        LOG_CONSOLE("\n Optimization successful!\n");
        LOG_CONSOLE("Estimated pose:\n");
        LOG_CONSOLE("Translation (mm): [%.7f, %.7f, %.7f]\n", t_est[0] * 1e3, t_est[1] * 1e3, t_est[2] * 1e3);
        LOG_CONSOLE("Rotation (mrad):  [%.7f, %.7f, %.7f]\n", xi_result[3] * 1e3, xi_result[4] * 1e3, xi_result[5] * 1e3);
        
        LOG_CONSOLE("\nTrue pose:\n");
        LOG_CONSOLE("Translation (mm): [%.7f, %.7f, %.7f]\n", tx * 1e3, ty * 1e3, tz * 1e3);
        LOG_CONSOLE("Rotation (mrad):  [%.7f, %.7f, %.7f]\n", theta_x * 1e3, theta_y * 1e3, theta_z * 1e3);
        
        // 计算误差
        lie_scalar_t trans_error[3] = {t_est[0] - tx, t_est[1] - ty, t_est[2] - tz};
        lie_scalar_t rot_error[3] = {xi_result[3] - theta_x, xi_result[4] - theta_y, xi_result[5] - theta_z};
        
        LOG_CONSOLE("\nErrors:\n");
        LOG_CONSOLE("Translation error (nm): [%.2e, %.2e, %.2e]\n", 
               trans_error[0]*1e9, trans_error[1]*1e9, trans_error[2]*1e9);
        LOG_CONSOLE("Rotation error (nrad):  [%.2e, %.2e, %.2e]\n", 
               rot_error[0]*1e9, rot_error[1]*1e9, rot_error[2]*1e9);

        LOG_CONSOLE("\nMaximum translation error (nm): [%.3e, %.3e, %.3e]\n", 
               max_trans_error[0]*1e9, max_trans_error[1]*1e9, max_trans_error[2]*1e9);
        LOG_CONSOLE("Maximum rotation error (nrad):  [%.3e, %.3e, %.3e]\n", 
               max_rot_error[0]*1e9, max_rot_error[1]*1e9, max_rot_error[2]*1e9);
        
        LOG_CONSOLE("\nFinal residuals (nm):   [%.2e, %.2e, %.2e]\n", 
                final_residual[0]*1e9, final_residual[1]*1e9, final_residual[2]*1e9);
        
        // 检查是否成功收敛到真实值
        lie_scalar_t trans_error_norm = sqrt(trans_error[0]*trans_error[0] + trans_error[1]*trans_error[1] + trans_error[2]*trans_error[2]);
        lie_scalar_t rot_error_norm = sqrt(rot_error[0]*rot_error[0] + rot_error[1]*rot_error[1] + rot_error[2]*rot_error[2]);
        
        if (trans_error_norm < 1e-9 && rot_error_norm < 1e-8) {
            LOG_CONSOLE("\n Perfect convergence achieved!\n");
        } else if (trans_error_norm < 1e-6 && rot_error_norm < 1e-6) {
            LOG_CONSOLE("\n Good convergence achieved!\n");
        } else {
            LOG_CONSOLE("\n  Convergence may be suboptimal. Consider adjusting parameters.\n");
        }
        
        LOG_CONSOLE("\n  Execution time: %.6f us\n", cpu_time_used);
        LOG_CONSOLE("Maximum execution time: %.6f us\n", max_exec_time);
    } else {
        LOG_CONSOLE(" Optimization failed with status: %d\n", status);
        LOG_CONSOLE("Execution time: %.6f us\n", cpu_time_used);
        LOG_CONSOLE("Status codes:\n");
        LOG_CONSOLE("  0 = Success\n");
        LOG_CONSOLE(" -1 = Invalid input\n");
        LOG_CONSOLE(" -2 = Memory allocation failed\n");
        LOG_CONSOLE(" -3 = Matrix inversion failed\n");
        LOG_CONSOLE(" -4 = Lambda too large\n");
        LOG_CONSOLE("\n  Execution time: %.6f us\n", cpu_time_used);
    }
    
    return 0;
}
#endif
