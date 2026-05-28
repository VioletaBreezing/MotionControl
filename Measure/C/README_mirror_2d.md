# 2D Mirror Reflection Simulation - C Implementation

## 概述

这是 `mirror_reflection_simulation_2d.m` MATLAB程序的C语言实现版本。该程序用于模拟和计算2D镜面反射系统中的光程距离，并使用梯度下降法优化3自由度位姿（x, y, rz）。

## 文件说明

- **mirror_reflection_simulation_2d.h**: 头文件，包含数据结构定义和函数声明
- **mirror_reflection_simulation_2d.c**: 主实现文件，包含所有核心算法
- **test_mirror_2d.c**: 测试程序，演示如何使用这些函数

## 主要功能

### 1. cal_light_distance
计算光程距离：h = f .* g

### 2. update_params
更新参数结构体中的中间变量，包括gamma、zeta、beta、sin_zeta、cos_zeta等。

### 3. cal_gradient
计算梯度向量，用于后续的优化过程。

### 4. cal_3dof
通过求解线性方程组计算3自由度位姿（x, y, rz）。

### 5. cal_3dof_grad
使用梯度下降法迭代优化3自由度位姿，包含8个步骤的时间统计。

### 6. mirror_sim_2d_main
主测试函数，运行N次随机测试并统计最大误差和计算时间。

## 编译方法

### Windows (使用MSVC)
```batch
cl /O2 test_mirror_2d.c mirror_reflection_simulation_2d.c /Fe:test_mirror_2d.exe
```

### Linux/Mac (使用GCC)
```bash
gcc -O2 -o test_mirror_2d test_mirror_2d.c mirror_reflection_simulation_2d.c -lm
```

### 启用浮点精度
如果需要单精度浮点数（float）而非双精度（double），在编译时添加宏定义：
```bash
gcc -O2 -DLIE_SOPHUS_USE_FLOAT -o test_mirror_2d test_mirror_2d.c mirror_reflection_simulation_2d.c -lm
```

### DSP平台优化 (C6678)
如果目标平台是TI C6678 DSP，添加相应的宏定义：
```bash
gcc -O2 -DSOC_C6678 -o test_mirror_2d test_mirror_2d.c mirror_reflection_simulation_2d.c -lm
```

## 使用方法

### 基本用法

```c
#include "mirror_reflection_simulation_2d.h"

int main() {
    // 运行10000次迭代测试
    mirror_sim_2d_main(10000);
    return 0;
}
```

### 单独使用各个函数

```c
#include "mirror_reflection_simulation_2d.h"

int main() {
    MirrorSim2D_Params params;
    
    // 初始化参数
    params.phi[0] = 1e-7; params.phi[1] = 2e-7; params.phi[2] = 3e-7;
    params.psi[0] = 1e-7; params.psi[1] = 2e-7; params.psi[2] = 3e-7;
    params.a_x[0] = -0.2; params.a_y[0] = 0.1;
    params.a_x[1] = -0.1; params.a_y[1] = -0.2;
    params.a_x[2] = 0.1;  params.a_y[2] = -0.2;
    params.b[0] = 0.1; params.b[1] = 0.1; params.b[2] = 0.1;
    
    lie_scalar_t x = 0.0, y = 0.0, rz = 1e-5;
    
    // 更新参数
    update_params(x, y, rz, &params);
    
    // 计算光程
    lie_scalar_t h[3];
    cal_light_distance(&params, h);
    
    printf("Light distances: h1=%.6e, h2=%.6e, h3=%.6e\n", h[0], h[1], h[2]);
    
    // 使用梯度下降法求解位姿
    lie_scalar_t x_est, y_est, rz_est, t[8];
    cal_3dof_grad(h, &params, rz, &x_est, &y_est, &rz_est, t);
    
    printf("Estimated pose: x=%.6e, y=%.6e, rz=%.6e\n", x_est, y_est, rz_est);
    
    return 0;
}
```

## 与MATLAB版本的对应关系

| MATLAB函数 | C函数 | 说明 |
|-----------|-------|------|
| cal_light_distance | cal_light_distance | 计算光程距离 |
| update_para | update_params | 更新参数 |
| cal_gradient | cal_gradient | 计算梯度 |
| cal_3dof | cal_3dof | 求解3自由度位姿 |
| cal_3dof_grad | cal_3dof_grad | 梯度下降优化 |
| main_func | mirror_sim_2d_main | 主测试函数 |

## 注意事项

1. **索引差异**: MATLAB使用1-based索引，C使用0-based索引。代码中已做相应调整。

2. **随机数生成**: C版本使用标准库的`rand()`函数，与MATLAB的`rand()`分布可能略有不同。

3. **时间测量**: C版本使用`clock()`函数进行时间测量，精度可能因平台而异。

4. **数值精度**: 默认使用double精度。如需float精度，定义`LIE_SOPHUS_USE_FLOAT`宏。

5. **除零保护**: 在计算f = beta/(2*beta^2-1)时，添加了分母接近零的保护。

## 性能特点

- 代码针对嵌入式系统和DSP平台进行了优化
- 支持SOC_C6678平台的代码段优化指令
- 最小化动态内存分配
- 使用内联数学运算提高性能

## 输出示例

运行测试程序后，将输出：
- 每个步骤的平均执行时间（微秒）
- 总平均计算时间
- 最大误差（x, y, rz方向）
- 最大计算时间

## 许可证

本项目遵循原MotionControl项目的许可协议。
