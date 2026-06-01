#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// 添加max宏定义
#ifndef max
#define max(a,b) ((a) > (b) ? (a) : (b))
#endif

#ifndef min
#define min(a,b) ((a) < (b) ? (a) : (b))
#endif

#ifdef _WIN32
#include <time.h>
#include <windows.h>
typedef LARGE_INTEGER timer_t;
timer_t freq;
#define timer_start(x)  do{QueryPerformanceCounter(&x);}while(0)
#define timer_stop(x)   do{QueryPerformanceCounter(&x);}while(0)
#define cal_time_us(x, y, t)  do{t = (double)(y.QuadPart - x.QuadPart) * 1000000.0 / freq.QuadPart;}while(0)
// #define random_seed (unsigned int)time(NULL)
#define random_seed (42)
#define my_sin  sin
#define my_cos  cos
// #define my_sin  custom_sin
// #define my_cos  custom_cos
#define my_atan atan
#define my_fabs fabs
#define my_sqrt sqrt
#define my_solver   linalg_solve_Axb_3x3
#elif defined(SOC_C6678)
#include "BSP_int.h"
#include "BSP_timer.h"
typedef double timer_t;
#define timer_start(x)  do{x = BSP_get_time();}while(0)
#define timer_stop(x)   do{x = BSP_get_time();}while(0)
#define cal_time_us(x, y, t) do{t = BSP_calc_us_interval(x, y);}while(0)
#define random_seed (42)
#define my_sin  sin
#define my_cos  cos
#define my_atan atan
#define my_fabs fabs
#define my_solver   linalg_solve_Axb_3x3
#endif

#ifdef SOC_C6678
#pragma CODE_SECTION(custom_sin, ".sa_code")
#pragma CODE_SECTION(custom_cos, ".sa_code")
#pragma CODE_SECTION(self_test, ".sa_code")
#pragma CODE_SECTION(sensor_model, ".sa_code")
#pragma CODE_SECTION(linalg_solve_Axb_3x3, ".sa_code")
#pragma CODE_SECTION(det3, ".sa_code")
#endif

static double PI = 3.14159265358979323846;

typedef struct {
    double phi[3];      // phi参数
    double psi[3];      // psi参数
    double a_x[3];      // a的x坐标
    double a_y[3];      // a的y坐标
    double b[3];        // b参数
} IFM_param_t;

typedef struct {
    double x, y, rz;
} IFM_state_t;

static double randn();
static double my_rms_func(double* arr, int N);
static void linalg_solve_Axb_3x3(double* A, double* b, double* x);
static inline double custom_sin(double x);
static inline double custom_cos(double x);
void sensor_model(IFM_state_t *pos, double *h, double *J, IFM_param_t *param)
{
    int i = 0;
    memset(h, 0, sizeof(double) * 3);
    memset(J, 0, sizeof(double) * 9);

    double zeta[3], gamma[3];

    for (i = 0; i < 3; i++) {
        zeta[i]  = pos->rz + param->psi[i];
        gamma[i] = pos->rz + param->psi[i] - param->phi[i];
    }
    // zeta[0] = zeta[0] - PI/2;

    for (i = 0; i < 1; i++) {
        double sin_gamma = my_sin(gamma[i]);
        double cos_gamma = my_cos(gamma[i]);
        double cos_2gamma = my_cos(2 * gamma[i]);
        double f = -cos_gamma / cos_2gamma;

        double sin_zeta = -my_cos(zeta[i]);
        double cos_zeta = my_sin(zeta[i]);

        double dx = pos->x - param->a_x[i];
        double dy = pos->y - param->a_y[i];

        double g = param->b[i] + dx * sin_zeta - dy * cos_zeta;
        h[i] = f * g;

        /* --- Jacobian --- */
        double df_drz = -sin_gamma * (cos_2gamma + 2) / (cos_2gamma * cos_2gamma);
        double dg_dx = sin_zeta;
        double dg_dy = -cos_zeta;
        double dg_drz = dx * cos_zeta + dy * sin_zeta;
        double dh_dx = f * dg_dx;
        double dh_dy = f * dg_dy;
        double dh_drz = df_drz * g + f * dg_drz;

        J[i*3] = dh_dx; J[i*3+1] = dh_dy; J[i*3+2] = dh_drz;
    }

    for (i = 1; i < 3; i++) {
        double sin_gamma = my_sin(gamma[i]);
        double cos_gamma = my_cos(gamma[i]);
        double cos_2gamma = my_cos(2 * gamma[i]);
        double f = -cos_gamma / cos_2gamma;

        double sin_zeta = my_sin(zeta[i]);
        double cos_zeta = my_cos(zeta[i]);

        double dx = pos->x - param->a_x[i];
        double dy = pos->y - param->a_y[i];

        double g = param->b[i] + dx * sin_zeta - dy * cos_zeta;
        h[i] = f * g;

        /* --- Jacobian --- */
        double df_drz = -sin_gamma * (cos_2gamma + 2) / (cos_2gamma * cos_2gamma);
        double dg_dx = sin_zeta;
        double dg_dy = -cos_zeta;
        double dg_drz = dx * cos_zeta + dy * sin_zeta;
        double dh_dx = f * dg_dx;
        double dh_dy = f * dg_dy;
        double dh_drz = df_drz * g + f * dg_drz;

        J[i*3] = dh_dx; J[i*3+1] = dh_dy; J[i*3+2] = dh_drz;
    }
}

int self_test(void) {
    int i = 0;

    /* --- 系统参数 ---*/
    static IFM_param_t param;
    param.a_x[0] = -386.5e-3; param.a_x[1] = -53.2e-3; param.a_x[2] = 53.8e-3;
    param.a_y[0] = 157.7e-3; param.a_y[1] = -241e-3; param.a_y[2] = -241e-3;

    param.b[0] = 213e-3; param.b[1] = 142e-3; param.b[2] = 142e-3;

    param.psi[0] = 1e-4; param.psi[1] = 1e-5; param.psi[2] = -1e-5;
    param.phi[0] = -2e-4; param.phi[1] = 2e-5; param.phi[2] = -2e-5;

    /* --- 性能评估 --- */
#ifdef _WIN32
    QueryPerformanceFrequency(&freq);
#endif
    timer_t start_tic, end_est_tic, end_precal_tic;
    double max_est_time = -1, max_precal_time = -1, avg_time, used_est_time, used_precal_time, total_time = 0;

    /* --- 轨迹设置 --- */
    double fs = 50;
    double Ts = 1 / fs;
    int N = (int)fs * 10;

    /* True trajectory (1 Hz sine) */
    double* x_true = (double*)malloc(sizeof(double) * N);
    double* y_true = (double*)malloc(sizeof(double) * N);
    double* theta_true = (double*)malloc(sizeof(double) * N);
    for (i = 0; i < N; i++) {
        x_true[i] = 2e-3 * sin(2 * PI * i * Ts);
        y_true[i] = 1e-1 + 1e-1 * sin(2 * PI * i * Ts);
        theta_true[i] = 2e-3 * sin(2 * PI * i * Ts);
    }

    /* Preallocate */
    double* x_est = (double*)malloc(sizeof(double) * N);
    double* y_est = (double*)malloc(sizeof(double) * N);
    double* theta_est = (double*)malloc(sizeof(double) * N);
    double* err_x = (double*)malloc(sizeof(double) * N);
    double* err_y = (double*)malloc(sizeof(double) * N);
    double* err_theta = (double*)malloc(sizeof(double) * N);
    double  h_true[3], h_meas[3], h_est[3], J[9];
    double x0, y0, theta0;

    /* --- Sensor noise (typical high-end interferometer) --- */
    double noise_rms = 0.31e-9 / 6;
    srand(random_seed);

    /* --- Initial guess --- */
    x_est[0] = 0.0;
    y_est[0] = 0.0;
    theta_est[0] = 0.0;

    printf("Running closed-loop estimation with %.2f nm RMS noise...\n", noise_rms*1e9);

    int k = 0;
    for (k = 0; k < N; k++) {
        /* Generate true measurement */
        IFM_state_t pos;
        pos.x = x_true[k]; pos.y = y_true[k]; pos.rz = theta_true[k];
        sensor_model(&pos, h_true, J, &param);

        /* Add realistic sensor noise */
        h_meas[0] = h_true[0] + noise_rms * randn();
        h_meas[1] = h_true[1] + noise_rms * randn();
        h_meas[2] = h_true[2] + noise_rms * randn();

        /* Initial guess: use previous estimate */
        if (k == 0) {
            x0 = x_est[0];
            y0 = y_est[0];
            theta0 = theta_est[0];
        } else {
            x0 = x_est[k-1];
            y0 = y_est[k-1];
            theta0 = theta_est[k-1];
        }

        /* --- 2-step Newton-Raphson --- */
        pos.x = x0; pos.y = y0; pos.rz = theta0;
        int iter = 0; double r[3], delta[3];

        timer_start(start_tic);
        for (iter = 0; iter < 2; iter++) {
            sensor_model(&pos, h_est, J, &param);
            r[0] = h_meas[0] - h_est[0];
            r[1] = h_meas[1] - h_est[1];
            r[2] = h_meas[2] - h_est[2];
            my_solver(J, r, delta);
            pos.x += delta[0];
            pos.y += delta[1];
            pos.rz += delta[2];
        }

        /* Record time cost */
        timer_stop(end_est_tic);
        cal_time_us(start_tic, end_est_tic, used_est_time);
        max_est_time = max(max_est_time, used_est_time);
        total_time += used_est_time;

        /* Store estimate */
        x_est[k] = pos.x;
        y_est[k] = pos.y;
        theta_est[k] = pos.rz;

        /* Compute error */
        err_x[k] = pos.x - x_true[k];
        err_y[k] = pos.y - y_true[k];
        err_theta[k] = pos.rz- theta_true[k];

        double err_x_debug = err_x[k];
        double err_y_debug = err_y[k];
        double err_theta_debug = err_theta[k];
    }

    /* --- Final analysis: use RMS error --- */
    double rms_err_x = my_rms_func(err_x, N);
    double rms_err_y = my_rms_func(err_y, N);
    double rms_err_theta = my_rms_func(err_theta, N);
    avg_time = total_time / N;

    double baseline = my_fabs(param.a_x[2] - param.a_x[1]);
    double equiv_noise_theta = noise_rms / baseline; // rad

    printf("\n✅ FINAL RESULT (with %.2f nm RMS sensor noise):\n", noise_rms*1e9);
    printf("RMS |Δx| = %.3f nm\n", rms_err_x*1e9);
    printf("RMS |Δy| = %.3f nm\n", rms_err_y*1e9);
    printf("Equivalent angle noise floor: %.2f nrad\n", equiv_noise_theta*1e9);
    printf("Achieved RMS Δθ             : %.2f nrad\n", rms_err_theta*1e9);
    printf("Maximum est time: %.3f us\n", max_est_time);
    printf("Average est time: %.3f us\n", avg_time);

    if (rms_err_x < 2*noise_rms && \
    rms_err_y < 2*noise_rms && \
    rms_err_theta < 2*equiv_noise_theta) {
        printf("\n🎉 SUCCESS: All errors near theoretical limit!\n");
    } else {
        printf("\n⚠️  Check model or noise assumptions.\n");
    }

    /* Free menmory */
    free(x_true);
    free(y_true);
    free(theta_true);

    free(x_est);
    free(y_est);
    free(theta_est);
    free(err_x);
    free(err_y);
    free(err_theta);
}

static double randn() {
    static int haveSpare = 0;
    static double spare;

    if (haveSpare) {
        haveSpare = 0;
        return spare;
    }

    haveSpare = 1;
    double u, v, s;
    do {
        u = (rand() / ((double) RAND_MAX)) * 2.0 - 1.0; // [-1, 1)
        v = (rand() / ((double) RAND_MAX)) * 2.0 - 1.0; // [-1, 1)
        s = u * u + v * v;
    } while (s >= 1.0 || s == 0.0);

    double mul = sqrt(-2.0 * log(s) / s);
    spare = v * mul;
    return u * mul;
}

static double my_rms_func(double* arr, int N) {
    double square_sum = 0.0;
    int i = 0;
    for (i = 0; i < N; i++) {
        square_sum += arr[i] * arr[i];
    }

    return sqrt(square_sum / N);
}

static inline double det3(const double* M) {
    // M is 3x3 row-major: [m00, m01, m02, m10, m11, m12, m20, m21, m22]
    return M[0] * (M[4] * M[8] - M[5] * M[7])
         - M[1] * (M[3] * M[8] - M[5] * M[6])
         + M[2] * (M[3] * M[7] - M[4] * M[6]);
}
static inline void linalg_solve_Axb_3x3(double* A, double* b, double* x) {
    double detA = det3(A);
    
    // 可选：检查 detA 是否为零（此处略去错误处理）
    // if (fabs(detA) < 1e-12) { /* handle singular matrix */ }

    // 构造 A_x: 第0列替换为 b
    double Ax[9] = {b[0], A[1], A[2],
                    b[1], A[4], A[5],
                    b[2], A[7], A[8]};
    // 构造 A_y: 第1列替换为 b
    double Ay[9] = {A[0], b[0], A[2],
                    A[3], b[1], A[5],
                    A[6], b[2], A[8]};
    // 构造 A_z: 第2列替换为 b
    double Az[9] = {A[0], A[1], b[0],
                    A[3], A[4], b[1],
                    A[6], A[7], b[2]};

    x[0] = det3(Ax) / detA;
    x[1] = det3(Ay) / detA;
    x[2] = det3(Az) / detA;
}

static inline double custom_sin(double x) {
    return x;
}

static inline double custom_cos(double x) {
    return 1;// - x * x / 2;
}

#ifdef _WIN32
void main() {
    self_test();

    double A[9] = {2, -1, 0,
                   -1, 2, -1,
                   0, -1, 2};
    double b[3] = {1, 0, 2};
    double x[3];
    linalg_solve_Axb_3x3(A, b, x);
    printf("x = [%f, %f, %f]\n", x[0], x[1], x[2]); // 应输出 [1, 1, 1]
}
#endif
