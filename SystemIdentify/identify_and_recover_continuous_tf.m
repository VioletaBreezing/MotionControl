clear; clc;

%% 1. 定义真实的连续系统（用于生成数据）
s = tf('s');
G_true = 10 / (s^2 + 2*s + 10);   % 二阶系统：ωn=√10, ζ=0.316
Ts = 0.001;                         % 采样周期（秒）

%% 2. 生成仿真数据
t = 0:Ts:10;                      % 时间向量
u = idinput(length(t), 'prbs', [], [0 1]);  % PRBS 输入信号（充分激励）
y_true = lsim(G_true, u, t);      % 理想输出
y = y_true + 0.02*randn(size(y_true));  % 添加小噪声

% 创建 iddata 对象（系统辨识所需格式）
data = iddata(y, u, Ts);

%% 3. 辨识离散 ARX 模型（或其他结构）
na = 2; nb = 2; nk = 1;           % 模型阶数：A(q)y = B(q)u
sysd = arx(data, [na nb nk]);     % 辨识离散模型

disp('辨识得到的离散模型：');
sysd

%% 4. 最大保真地转换为连续传递函数（关键步骤！）
% 先转为状态空间（更稳定）
sysd_ss = ss(sysd);
% 使用 ZOH 逆变换（假设辨识时输入由零阶保持器生成）
sysc_ss = d2c(sysd_ss, 'zoh');    
% 转为传递函数形式
sysc_tf = tf(sysc_ss);

disp('还原的连续传递函数：');
sysc_tf

%% 5. 验证保真度（修正版）
% 将还原的连续模型重新离散化
sysd_check = c2d(sysc_ss, Ts, 'zoh');

% 计算系统差异（H-infinity 范数）
diff_norm = norm(sysd - sysd_check, 2);
fprintf('辨识模型与回验模型的 H∞ 范数误差: %.2e\n', diff_norm);

% 时域验证：用相同输入仿真
[y_id, ~] = lsim(sysd, u, t);
[y_check, ~] = lsim(sysd_check, u, t);
mse_check = mean((y_id - y_check).^2);
fprintf('离散模型回验 MSE: %.2e\n', mse_check);

% 绘图：频率响应
figure;
bode(G_true, 'k--', sysc_tf, 'b-', {0.1, 100});
legend('真实系统', '辨识+还原系统');
title('频率响应对比（验证保真度）');
grid on;

% 时域响应
figure;
plot(t, y_true, 'k--', t, y, 'r.', t, lsim(sysc_tf, u, t), 'b-', 'LineWidth', 1.2);
legend('真实输出', '含噪测量', '还原连续模型预测');
xlabel('时间 (s)'); ylabel('输出 y(t)');
title('时域响应对比');
grid on;