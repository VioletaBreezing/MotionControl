% 定义符号变量
syms a11 a12 a21 a22 a31 a32 real

% 构造矩阵 A (3x2)
A = [a11, a12;
     a21, a22;
     a31, a32];

% 计算 (A' * A)^(-1) * A'
pseudo_inv = inv(A' * A) * A';

% 显示结果
disp('Symbolic expression of (A^T A)^{-1} A^T:');
pretty(pseudo_inv)

% 可选：简化表达式
pseudo_inv_simplified = simplify(pseudo_inv);
disp('Simplified version:');
pretty(pseudo_inv_simplified)