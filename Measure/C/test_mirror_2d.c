#include "mirror_reflection_simulation_2d.h"
#include <stdio.h>

/**
 * @brief 主函数 - 测试2D镜面反射仿真
 */
int main(void) {
    printf("========================================\n");
    printf("Mirror Reflection Simulation 2D Test\n");
    printf("========================================\n\n");
    
    // 运行10000次迭代测试（与MATLAB代码一致）
    int N = 500000;
    mirror_sim_2d_main(N);
    
    printf("\n========================================\n");
    printf("Test completed successfully!\n");
    printf("========================================\n");
    
    return 0;
}
