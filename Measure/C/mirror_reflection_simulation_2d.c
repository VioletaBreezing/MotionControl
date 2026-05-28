#include "mirror_reflection_simulation_2d.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#ifdef _WIN32
#include <windows.h>
#endif

#ifdef SOC_C6678
#pragma CODE_SECTION(cal_light_distance, ".sa_code")
#pragma CODE_SECTION(update_params, ".sa_code")
#pragma CODE_SECTION(cal_gradient, ".sa_code")
#pragma CODE_SECTION(cal_3dof, ".sa_code")
#pragma CODE_SECTION(cal_3dof_grad, ".sa_code")
#pragma CODE_SECTION(mirror_sim_2d_main, ".sa_code")
#endif

static const lie_scalar_t PI = 3.14159265358979323846;

static inline lie_scalar_t my_sin(lie_scalar_t x)
{
    return x;// - x * x * x / 6; // + x * x * x * x * x / 120;
};

static inline lie_scalar_t my_cos(lie_scalar_t x)
{
    return 1;// - x * x / 2; // + x * x * x * x / 24;
};

static inline lie_scalar_t my_atan(lie_scalar_t x)
{ 
    return x;// - x * x * x / 3; // + x * x * x * x * x / 15;
}

/**
 * @brief 计算光程距离 h = f .* g
 */
static inline void cal_light_distance(const MirrorSim2D_Params* params, lie_scalar_t* h) {
    h[0] = params->f[0] * params->g[0];
    h[1] = params->f[1] * params->g[1];
    h[2] = params->f[2] * params->g[2];
}

/**
 * @brief 更新参数结构体中的中间变量
 */
static inline void update_params(lie_scalar_t x, lie_scalar_t y, lie_scalar_t rz, MirrorSim2D_Params* params) {
    int i;
    
    for (i = 0; i < 3; i++) {
        // gamma(i) = rz + psi(i) - phi(i)
        params->gamma[i] = rz + params->psi[i] - params->phi[i];
        
        // zeta(i) = rz + psi(i)
        params->zeta[i] = rz + params->psi[i];
        
        // // 对于i==0，需要减去pi/2（MATLAB索引从1开始，C从0开始）
        // if (i == 0) {
        //     params->zeta[i] = params->zeta[i] - PI / 2.0;
        // }
        if (i == 0)
        {
            params->sin_zeta[0] = -lie_cos(params->zeta[0]);
            params->cos_zeta[0] = lie_sin(params->zeta[0]);
        }
        else
        {
            // 计算sin和cos
            params->sin_zeta[i] = lie_sin(params->zeta[i]);
            params->cos_zeta[i] = lie_cos(params->zeta[i]);
        }

        // beta(i) = -cos(gamma(i))
        params->beta[i] = -lie_cos(params->gamma[i]);
        
        // f(i) = beta(i) / (2*beta(i)^2 - 1)
        // lie_scalar_t beta_sq = params->beta[i] * params->beta[i];
        // lie_scalar_t denom = 2.0 * beta_sq - 1.0;
        // if (lie_fabs(denom) < 1e-12) {
        //     params->f[i] = 0.0; // 避免除零
        // } else {
        //     params->f[i] = params->beta[i] / denom;
        // }
        params->f[i] = -1;
        
        // dx(i) = x - a{i}(1)
        params->dx[i] = x - params->a_x[i];
        
        // dy(i) = y - a{i}(2)
        params->dy[i] = y - params->a_y[i];
        
        // g(i) = b(i) + dx(i)*sin_zeta(i) - dy(i)*cos_zeta(i)
        params->g[i] = params->b[i] + params->dx[i] * params->sin_zeta[i] 
                      - params->dy[i] * params->cos_zeta[i];
    }
}

/**
 * @brief 计算梯度
 */
static inline void cal_gradient(const MirrorSim2D_Params* params, lie_scalar_t* grad) {
    int i;
    
    for (i = 0; i < 3; i++) {
        // lie_scalar_t beta = params->beta[i];
        // lie_scalar_t f = params->f[i];
        // lie_scalar_t sin_zeta = params->sin_zeta[i];
        // lie_scalar_t cos_zeta = params->cos_zeta[i];
        // lie_scalar_t g = params->g[i];
        // lie_scalar_t dx = params->dx[i];
        // lie_scalar_t dy = params->dy[i];
        
        // dfdr = -((2*beta^2+1)/(2*beta^2-1))^2 * sin(gamma(i))
        // lie_scalar_t beta_sq = beta * beta;
        // lie_scalar_t numer = 2.0 * beta_sq + 1.0;
        // lie_scalar_t denom = 2.0 * beta_sq - 1.0;
        
        // lie_scalar_t ratio;
        // if (lie_fabs(denom) < 1e-12) {
        //     ratio = 0.0;
        // } else {
        //     ratio = numer / denom;
        // }
        // lie_scalar_t dfdr = -(ratio * ratio) * lie_sin(params->gamma[i]);
        lie_scalar_t dfdr = -9.0 * params->gamma[i];
        
        // dgdr = dx * cos_zeta + dy * sin_zeta
        // lie_scalar_t dgdr = dx * cos_zeta + dy * sin_zeta;
        lie_scalar_t dgdr = params->dx[i] * params->cos_zeta[i] + params->dy[i] * params->sin_zeta[i];
        
        // grad(i) = dfdr * g + f * dgdr
        // grad[i] = dfdr * g + f * dgdr;
        grad[i] = dfdr * params->g[i] - dgdr;
    }
}

/**
 * @brief 计算3自由度位姿 (x, y, rz)
 */
static inline void cal_3dof(const lie_scalar_t* h_measure, const MirrorSim2D_Params* params, 
              lie_scalar_t rz0, int flag, 
              lie_scalar_t* tx, lie_scalar_t* ty, lie_scalar_t* rz_out) {
    lie_scalar_t Ly = params->a_x[1] - params->a_x[0]; // a{2}(1) - a{1}(1)
    lie_scalar_t rz;
    MirrorSim2D_Params temp_params;
    
    // 复制参数以避免修改原始数据
    memcpy(&temp_params, params, sizeof(MirrorSim2D_Params));
    
    if (flag == 1) {
        rz = rz0;
    } else {
        // rz = atan((h_measure(3)-h_measure(2))/Ly)
        lie_scalar_t atg = (h_measure[2] - h_measure[1]) / Ly;
        rz = lie_atan(atg);
        // 更新参数
        update_params(0.0, 0.0, rz, &temp_params);
    }
    
    // 构建线性方程组 U * d = v
    // U是3x2矩阵，v是3x1向量
    lie_scalar_t U[6]; // 3x2矩阵，按行存储
    lie_scalar_t v[3];
    
    int i;
    for (i = 0; i < 3; i++) {
        lie_scalar_t b_val = temp_params.b[i];
        lie_scalar_t ax = temp_params.a_x[i];
        lie_scalar_t ay = temp_params.a_y[i];
        lie_scalar_t f_val = temp_params.f[i];
        lie_scalar_t sin_z = temp_params.sin_zeta[i];
        lie_scalar_t cos_z = temp_params.cos_zeta[i];
        
        // U(i,1) = f * sin_zeta; U(i,2) = -f * cos_zeta
        U[i * 2 + 0] = f_val * sin_z;
        U[i * 2 + 1] = -f_val * cos_z;
        
        // v(i) = f * (b - ax*sin_zeta + ay*cos_zeta)
        v[i] = f_val * (b_val - ax * sin_z + ay * cos_z);
    }
    
    // h_hat = h_measure - v
    lie_scalar_t h_hat[3];
    for (i = 0; i < 3; i++) {
        h_hat[i] = h_measure[i] - v[i];
    }
    
    // 求解最小二乘问题: d = U \ h_hat
    // 使用正规方程: (U'*U) * d = U' * h_hat
    // U'*U 是 2x2 矩阵
    lie_scalar_t UtU[4]; // 2x2矩阵
    UtU[0] = U[0]*U[0] + U[2]*U[2] + U[4]*U[4]; // sum(U(:,1).^2)
    UtU[1] = U[0]*U[1] + U[2]*U[3] + U[4]*U[5]; // sum(U(:,1).*U(:,2))
    UtU[2] = UtU[1];                             // 对称
    UtU[3] = U[1]*U[1] + U[3]*U[3] + U[5]*U[5]; // sum(U(:,2).^2)
    
    // U' * h_hat
    lie_scalar_t Uth[2];
    Uth[0] = U[0]*h_hat[0] + U[2]*h_hat[1] + U[4]*h_hat[2];
    Uth[1] = U[1]*h_hat[0] + U[3]*h_hat[1] + U[5]*h_hat[2];
    
    // 求解 2x2 线性系统
    lie_scalar_t det = UtU[0] * UtU[3] - UtU[1] * UtU[2];
    if (lie_fabs(det) < 1e-12) {
        *tx = 0.0;
        *ty = 0.0;
    } else {
        *tx = (UtU[3] * Uth[0] - UtU[1] * Uth[1]) / det;
        *ty = (-UtU[2] * Uth[0] + UtU[0] * Uth[1]) / det;
    }
    
    *rz_out = rz;
}

/**
 * @brief 使用梯度下降法计算3自由度位姿
 */
static inline void cal_3dof_grad(const lie_scalar_t* h, MirrorSim2D_Params* params, lie_scalar_t rz0,
                   lie_scalar_t* tx, lie_scalar_t* ty, lie_scalar_t* rz, lie_scalar_t* t) {
    // 初始化时间统计
    int i;
    for (i = 0; i < 8; i++) {
        t[i] = 0.0;
    }
    
#ifdef _WIN32
    LARGE_INTEGER start, end, freq;
    QueryPerformanceFrequency(&freq);  // 获取计数器频率
#else
    clock_t start, end;
#endif
    lie_scalar_t x, y;
    
    // 第1步：cal_3dof初始估计
#ifdef _WIN32
    QueryPerformanceCounter(&start);
#else
    start = clock();
#endif
    cal_3dof(h, params, rz0, 0, &x, &y, &rz0);
#ifdef _WIN32
    QueryPerformanceCounter(&end);
    t[0] += (lie_scalar_t)(end.QuadPart - start.QuadPart) * 1000000.0 / freq.QuadPart;  // 微秒
#else
    end = clock();
    t[0] += (lie_scalar_t)(end - start) / CLOCKS_PER_SEC;
#endif
    
    // 迭代优化（只迭代1次，与MATLAB代码一致）
    for (int iter = 0; iter < 3; iter++) {
        // 第2步：update_params
#ifdef _WIN32
        QueryPerformanceCounter(&start);
#else
        start = clock();
#endif
        update_params(x, y, rz0, params);
#ifdef _WIN32
        QueryPerformanceCounter(&end);
        t[1] += (lie_scalar_t)(end.QuadPart - start.QuadPart) * 1000000.0 / freq.QuadPart;  // 微秒
#else
        end = clock();
        t[1] += (lie_scalar_t)(end - start) / CLOCKS_PER_SEC;
#endif
        
        // 第3步：cal_gradient
#ifdef _WIN32
        QueryPerformanceCounter(&start);
#else
        start = clock();
#endif
        lie_scalar_t grad[3];
        cal_gradient(params, grad);
#ifdef _WIN32
        QueryPerformanceCounter(&end);
        t[2] += (lie_scalar_t)(end.QuadPart - start.QuadPart) * 1000000.0 / freq.QuadPart;  // 微秒
#else
        end = clock();
        t[2] += (lie_scalar_t)(end - start) / CLOCKS_PER_SEC;
#endif
        
        // 第4步：计算残差 residual = cal_light_distance(params) - h
#ifdef _WIN32
        QueryPerformanceCounter(&start);
#else
        start = clock();
#endif
        lie_scalar_t h_pred[3];
        cal_light_distance(params, h_pred);
        lie_scalar_t residual[3];
        // for (i = 0; i < 3; i++) {
        //     residual[i] = h_pred[i] - h[i];
        // }
        residual[0] = h_pred[0] - h[0];
        residual[1] = h_pred[1] - h[1];
        residual[2] = h_pred[2] - h[2];
#ifdef _WIN32
        QueryPerformanceCounter(&end);
        t[3] += (lie_scalar_t)(end.QuadPart - start.QuadPart) * 1000000.0 / freq.QuadPart;  // 微秒
#else
        end = clock();
        t[3] += (lie_scalar_t)(end - start) / CLOCKS_PER_SEC;
#endif
        
        // 第5步：计算更新量 dr = -50 * grad' * residual
#ifdef _WIN32
        QueryPerformanceCounter(&start);
#else
        start = clock();
#endif
        lie_scalar_t dr = -50.0 * (grad[0] * residual[0] + grad[1] * residual[1] + grad[2] * residual[2]);
#ifdef _WIN32
        QueryPerformanceCounter(&end);
        t[4] += (lie_scalar_t)(end.QuadPart - start.QuadPart) * 1000000.0 / freq.QuadPart;  // 微秒
#else
        end = clock();
        t[4] += (lie_scalar_t)(end - start) / CLOCKS_PER_SEC;
#endif
        
        // 第6步：更新rz
#ifdef _WIN32
        QueryPerformanceCounter(&start);
#else
        start = clock();
#endif
        rz0 = rz0 + dr;
#ifdef _WIN32
        QueryPerformanceCounter(&end);
        t[5] += (lie_scalar_t)(end.QuadPart - start.QuadPart) * 1000000.0 / freq.QuadPart;  // 微秒
#else
        end = clock();
        t[5] += (lie_scalar_t)(end - start) / CLOCKS_PER_SEC;
#endif
        
        // 第7步：再次update_params
#ifdef _WIN32
        QueryPerformanceCounter(&start);
#else
        start = clock();
#endif
        update_params(x, y, rz0, params);
#ifdef _WIN32
        QueryPerformanceCounter(&end);
        t[6] += (lie_scalar_t)(end.QuadPart - start.QuadPart) * 1000000.0 / freq.QuadPart;  // 微秒
#else
        end = clock();
        t[6] += (lie_scalar_t)(end - start) / CLOCKS_PER_SEC;
#endif
        
        // 第8步：再次cal_3dof
#ifdef _WIN32
        QueryPerformanceCounter(&start);
#else
        start = clock();
#endif
        cal_3dof(h, params, rz0, 1, &x, &y, &rz0);
#ifdef _WIN32
        QueryPerformanceCounter(&end);
        t[7] += (lie_scalar_t)(end.QuadPart - start.QuadPart) * 1000000.0 / freq.QuadPart;  // 微秒
#else
        end = clock();
        t[7] += (lie_scalar_t)(end - start) / CLOCKS_PER_SEC;
#endif
    }
    
    *tx = x;
    *ty = y;
    *rz = rz0;
}

/**
 * @brief 运行主测试函数
 */
void mirror_sim_2d_main(int N) {
    lie_scalar_t max_errx = -1.0;
    lie_scalar_t max_erry = -1.0;
    lie_scalar_t max_errrz = -1.0;
    lie_scalar_t max_time = -1.0;
    
#ifdef _WIN32
    LARGE_INTEGER total_start, total_end, freq;
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&total_start);
    lie_scalar_t t_overhead = (lie_scalar_t)(clock() - total_start.QuadPart) * 1000000.0 / freq.QuadPart;
#else
    clock_t total_start = clock();
    lie_scalar_t t_overhead = (lie_scalar_t)(clock() - total_start) / CLOCKS_PER_SEC;
#endif
    
    lie_scalar_t t_slice[8];
    lie_scalar_t t_slice_max[8];
    lie_scalar_t t_total = 0;
    memset(t_slice, 0, sizeof(t_slice));
    memset(t_slice_max, 0, sizeof(t_slice_max));
    
    // 初始化随机数种子
    srand((unsigned int)time(NULL));

    int i;
    lie_scalar_t *px, *py, *prz, *t_axis;
    px = (lie_scalar_t *)malloc(sizeof(lie_scalar_t) * N);
    py = (lie_scalar_t *)malloc(sizeof(lie_scalar_t) * N);
    prz = (lie_scalar_t *)malloc(sizeof(lie_scalar_t) * N);
    t_axis = (lie_scalar_t *)malloc(sizeof(lie_scalar_t) * N);
    lie_scalar_t Ts = 2e-4; // 200us采样时间
    for (i = 0; i < N; i++) { 
        px[i] = 1e-1 * sin(2*PI*i*Ts);
        py[i] = 2e-1 * cos(2*PI*i*Ts);
        prz[i] = 1e-3 * sin(2*PI*i*Ts);
        t_axis[i] = i * Ts;
    }
    
    int iter;
    for (iter = 0; iter < N; iter++) {
        MirrorSim2D_Params params;
        
        // 初始化参数
        // phi = [rand()*1e-6, rand()*1e-6, rand()*1e-6]
        params.phi[0] = -2e-3;// * 2.0 * (((lie_scalar_t)rand() / RAND_MAX) - 0.5);
        params.phi[1] = -2e-3;// * 2.0 * (((lie_scalar_t)rand() / RAND_MAX) - 0.5);
        params.phi[2] = -2e-3;// * 2.0 * (((lie_scalar_t)rand() / RAND_MAX) - 0.5);
        
        // psi = [rand()*1e-6, rand()*1e-6, rand()*1e-6]
        params.psi[0] = 2e-3;// * 2.0 * (((lie_scalar_t)rand() / RAND_MAX) - 0.5);
        params.psi[1] = 2e-3;// * 2.0 * (((lie_scalar_t)rand() / RAND_MAX) - 0.5);
        params.psi[2] = 2e-3;// * 2.0 * (((lie_scalar_t)rand() / RAND_MAX) - 0.5);
        
        // a = {[-0.2; 0.1], [-0.1; -0.2], [0.1; -0.2]}
        params.a_x[0] = -0.25; params.a_y[0] = 0.1;
        params.a_x[1] = -0.10; params.a_y[1] = -0.2;
        params.a_x[2] = 0.10;  params.a_y[2] = -0.2;
        
        // b = [0.1; 0.1; 0.1]
        params.b[0] = 0.12;
        params.b[1] = 0.125;
        params.b[2] = 0.124;
        
        // x = 0 * 2 * (rand()-0.5)
        // lie_scalar_t x = 1e-3 * 2.0 * (((lie_scalar_t)rand() / RAND_MAX) - 0.5);
        
        // // y = 0 * rand()
        // lie_scalar_t y = 2e-1 * ((lie_scalar_t)rand() / RAND_MAX);
        
        // // rz = 1e-4 * (rand()-0.5)
        // lie_scalar_t rz = 1e-3 * 2.0 * (((lie_scalar_t)rand() / RAND_MAX) - 0.5);
        lie_scalar_t x = px[iter];
        lie_scalar_t y = py[iter];
        lie_scalar_t rz = prz[iter];
        
        // 更新参数
        update_params(x, y, rz, &params);
        
        // 计算光程
        lie_scalar_t h[3];
        cal_light_distance(&params, h);
        
        // 使用梯度下降法求解
#ifdef _WIN32
        LARGE_INTEGER step_start, step_end;
        QueryPerformanceCounter(&step_start);
#else
        clock_t step_start = clock();
#endif
        lie_scalar_t xg, yg, rzg;
        lie_scalar_t t[8];
        cal_3dof_grad(h, &params, rz, &xg, &yg, &rzg, t);
#ifdef _WIN32
        QueryPerformanceCounter(&step_end);
        lie_scalar_t elapsed_time = (lie_scalar_t)(step_end.QuadPart - step_start.QuadPart) * 1000000.0 / freq.QuadPart;  // 微秒
#else
        clock_t step_end = clock();
        lie_scalar_t elapsed_time = (lie_scalar_t)(step_end - step_start) / CLOCKS_PER_SEC;
#endif
        
        // 累加时间统计
        int i;
        lie_scalar_t sum_t_slice = 0;
        for (i = 0; i < 8; i++) {
            sum_t_slice += t[i];
            if (t[i] > t_slice_max[i]) t_slice_max[i] = t[i];
        }
        if (sum_t_slice > t_total) 
        {
             for (i = 0; i < 8; i++) {
                t_slice[i] = t[i];
            }
            t_total = sum_t_slice;
        }
        
        // 计算误差
        lie_scalar_t errxg = xg - x;
        lie_scalar_t erryg = yg - y;
        lie_scalar_t errrzg = rzg - rz;
        
        // 更新最大误差
        lie_scalar_t abs_errxg = lie_fabs(errxg);
        lie_scalar_t abs_erryg = lie_fabs(erryg);
        lie_scalar_t abs_errrzg = lie_fabs(errrzg);
        
        if (abs_errxg > max_errx) max_errx = abs_errxg;
        if (abs_erryg > max_erry) max_erry = abs_erryg;
        if (abs_errrzg > max_errrz) max_errrz = abs_errrzg;
        if (elapsed_time > max_time) max_time = elapsed_time;
    }

    lie_scalar_t avg_time = max_time;  // max_time已经是微秒
    
    // 输出结果
    printf("Average time per step (microseconds):\n");
    printf("  Step 1 (cal_3dof):           %.3f us\n", t_slice[0]);
    printf("  Step 2 (update_params):      %.3f us\n", t_slice[1]);
    printf("  Step 3 (cal_gradient):       %.3f us\n", t_slice[2]);
    printf("  Step 4 (cal_residual):       %.3f us\n", t_slice[3]);
    printf("  Step 5 (calc_dr):            %.3f us\n", t_slice[4]);
    printf("  Step 6 (update_rz):          %.3f us\n", t_slice[5]);
    printf("  Step 7 (update_params):      %.3f us\n", t_slice[6]);
    printf("  Step 8 (cal_3dof):           %.3f us\n", t_slice[7]);
    printf("\nMaximum errors:\n");
    printf("  Max error x:                 %.6e\n", max_errx);
    printf("  Max error y:                 %.6e\n", max_erry);
    printf("  Max error rz:                %.6e\n", max_errrz);
    printf("  Max computation time:        %.3f us\n", max_time);

    printf("\nMaximum time per step (microseconds):\n");
    printf("  Step 1 (cal_3dof):           %.3f us\n", t_slice_max[0]);
    printf("  Step 2 (update_params):      %.3f us\n", t_slice_max[1]);
    printf("  Step 3 (cal_gradient):       %.3f us\n", t_slice_max[2]);
    printf("  Step 4 (cal_residual):       %.3f us\n", t_slice_max[3]);
    printf("  Step 5 (calc_dr):            %.3f us\n", t_slice_max[4]);
    printf("  Step 6 (update_rz):          %.3f us\n", t_slice_max[5]);
    printf("  Step 7 (update_params):      %.3f us\n", t_slice_max[6]);
    printf("  Step 8 (cal_3dof):           %.3f us\n", t_slice_max[7]);

    free(t_axis);
    free(px);
    free(py);
    free(prz);
}
