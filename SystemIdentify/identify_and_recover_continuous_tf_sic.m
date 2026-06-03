clear; clc;

%% 1. 定义真实的连续系统（用于生成数据）
s = tf('s');
G_true = 10 / (s^2 + 2*s + 10);   % 二阶系统：ωn=√10, ζ=0.316
Ts = 1e-4;                       % 采样周期（秒）

%% 2. 生成仿真数据
t = 0:Ts:10;                      % 时间向量（10秒）
u = idinput(length(t), 'prbs', [], [0 1]);  % PRBS 输入信号
y_true = lsim(G_true, u, t);      % 理想输出
y = y_true + 1e-6*randn(size(y_true));  % 添加小噪声

% 创建 iddata 对象（系统辨识所需格式）
data = iddata(y, u, Ts);

%% 3. 使用 arxstruc + selstruc 自动选择 ARX 模型阶数（基于 BIC）
% 定义搜索范围（注意：na >= 1, nb >= 1, nk >= 1）
na_range = 1:4;    % A 多项式阶数
nb_range = 1:4;    % B 多项式阶数
nk_range = 1;      % 输入延迟（通常设为 1，或 1:2）

% 生成所有可能的 [na nb nk] 组合
nn = struc(na_range, nb_range, nk_range);

% 使用 arxstruc 计算每个结构的拟合损失（使用同一数据进行估计和验证）
% 注意：第二个 data 是"验证集"，这里用训练集代替（也可分割数据）
V = arxstruc(data, data, nn);

% 方法 1：使用 selstruc 自动选择（支持 'aic', 'bic', 'mdl' 等）
best_nn_bic = selstruc(V, 'bic');   % 基于 BIC
best_nn_aic = selstruc(V, 'aic');   % 基于 AIC（可选）

fprintf('BIC 选择的模型阶数: na=%d, nb=%d, nk=%d\n', best_nn_bic(1), best_nn_bic(2), best_nn_bic(3));
fprintf('AIC 选择的模型阶数: na=%d, nb=%d, nk=%d\n', best_nn_aic(1), best_nn_aic(2), best_nn_aic(3));

% 选择 BIC 结果（更倾向于简约模型）
opt_na = best_nn_bic(1);
opt_nb = best_nn_bic(2);
opt_nk = best_nn_bic(3);

%% 4. 使用最优阶数辨识 ARX 模型
sysd = arx(data, [opt_na, opt_nb, opt_nk]);

disp('辨识得到的离散 ARX 模型（BIC 选阶）：');
sysd

%% 5. 转换为连续传递函数（关键步骤）
sysd_ss = ss(sysd);                % 转为状态空间（更稳定）
sysc_ss = d2c(sysd_ss, 'zoh');     % 零阶保持器逆变换
sysc_tf = tf(sysc_ss);             % 转为传递函数

disp('还原的连续传递函数：');
sysc_tf

%% 6. 验证保真度
% 将还原的连续模型重新离散化
sysd_check = c2d(sysc_ss, Ts, 'zoh');

% 计算 H-infinity 范数误差（或使用 2-范数）
diff_norm = norm(sysd - sysd_check, inf);
fprintf('辨识模型与回验模型的 H∞ 范数误差: %.2e\n', diff_norm);

% 时域回验 MSE
[y_id, ~] = lsim(sysd, u, t);
[y_check, ~] = lsim(sysd_check, u, t);
mse_check = mean((y_id - y_check).^2);
fprintf('离散模型回验 MSE: %.2e\n', mse_check);

%% 7. 绘图：频率响应对比
figure;
bode(G_true, 'k--', sysc_tf, 'b-', {0.1, 100});
legend('真实系统', '辨识+还原系统', 'Location', 'best');
title('频率响应对比（BIC 自动选阶）');
grid on;

% 时域响应对比
figure;
plot(t, y_true, 'k--', 'LineWidth', 1.2); hold on;
plot(t, y, 'r.', 'MarkerSize', 8);
plot(t, lsim(sysc_tf, u, t), 'b-', 'LineWidth', 1.2);
legend('真实输出', '含噪测量', '还原连续模型预测', 'Location', 'best');
xlabel('时间 (s)'); ylabel('输出 y(t)');
title('时域响应对比（BIC 自动选阶）');
grid on;

%% 可选：可视化结构选择结果（BIC/AIC 曲线）
% 显示损失 vs 参数数量
figure;
selstruc(V, 'bic');  % 弹出交互式选阶图（关闭即可）