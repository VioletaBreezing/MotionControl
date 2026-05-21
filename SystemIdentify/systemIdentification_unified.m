clear; close all; clc;
Ts = 4e-4;  % 控制周期
Fs = 1 / Ts;
fmin = 1; fmax = min(Fs/4, 400);

%% 维度设置
input_dim = 3;      % 输入维度
output_dim = 3;     % 输出维度
actuator_dim = 3;   % 作动器维度

input_dim_name = {"X", "Y", "RZ"};      % 各维度名称
output_dim_name = {"X", "Y", "RZ"};
actuator_dim_name = {"X", "Y1", "Y2"};

%% 数据导入
data_dir = "RS_SS_系统辨识_400u_AC0.29467/";
x_dir = "X_inject/";
y_dir = "Y_inject/";
rz_dir = "Rz_inject/";
pos_reverse_phase = true;

u1_file         = data_dir + x_dir + "20260512_170405_WSSS_CTRL_OUT_X.ini";
y1_u1_file      = data_dir + x_dir + "20260512_170405_WSSS_ERR_POS_X.ini";
y2_u1_file      = data_dir + x_dir + "20260512_170405_WSSS_ERR_POS_Y.ini";
y3_u1_file      = data_dir + x_dir + "20260512_170405_WSSS_ERR_POS_RZ.ini";
inject_u1_file  = data_dir + x_dir + '20260512_170405_WSSS_TEST_SIGNAL_X.ini';

u2_file         = data_dir + y_dir + "20260512_170828_WSSS_CTRL_OUT_Y.ini";
y1_u2_file      = data_dir + y_dir + "20260512_170828_WSSS_ERR_POS_X.ini";
y2_u2_file      = data_dir + y_dir + "20260512_170828_WSSS_ERR_POS_Y.ini";
y3_u2_file      = data_dir + y_dir + "20260512_170828_WSSS_ERR_POS_RZ.ini";
inject_u2_file  = data_dir + y_dir + '20260512_170828_WSSS_TEST_SIGNAL_Y.ini';

u3_file         = data_dir + rz_dir + "20260512_171455_WSSS_CTRL_OUT_RZ.ini";
y1_u3_file      = data_dir + rz_dir + "20260512_171455_WSSS_ERR_POS_X.ini";
y2_u3_file      = data_dir + rz_dir + "20260512_171455_WSSS_ERR_POS_Y.ini";
y3_u3_file      = data_dir + rz_dir + "20260512_171455_WSSS_ERR_POS_RZ.ini";
inject_u3_file  = data_dir + rz_dir + '20260512_171455_WSSS_TEST_SIGNAL_RZ.ini';

u_files = {u1_file, u2_file, u3_file}';
y_files = {y1_u1_file, y2_u1_file, y3_u1_file;
           y1_u2_file, y2_u2_file, y3_u2_file;
           y1_u3_file, y2_u3_file, y3_u3_file};
inject_files = {inject_u1_file, inject_u2_file, inject_u3_file}';

% 读文件
u_time_domain       = cell(input_dim, 1);
y_time_domain       = cell(input_dim, output_dim);
inject_time_domain  = cell(input_dim, 1);

for i = 1:input_dim
    u_time_domain{i}      = readSingleColumnData(u_files{i});
    inject_time_domain{i} = readSingleColumnData(inject_files{i});

    for j = 1:output_dim
        if pos_reverse_phase
            y_time_domain{i, j} = -readSingleColumnData(y_files{i, j});
        else
            y_time_domain{i, j} = readSingleColumnData(y_files{i, j});
        end
    end
end

%% 画时域图
for i = 1:input_dim
    displayTimedomain(u_time_domain{i}, y_time_domain{i, i}, ...
        "Title1", input_dim_name{i}  + "输入", ...
        "Title2", output_dim_name{i} + "输出", ...
        "YLabel1", "力/N", ...
        "YLabel2", "位移/m", ...
        "Ts", Ts);
end

%% 计算MIMO传递函数
sys.tf = cell(input_dim, output_dim);
sys.siso_plant_frd_cont = cell(input_dim, 1);
sys.siso_plant_frd_disc = cell(input_dim, 1);
sys.coh = cell(input_dim, output_dim);
sys.Fs = Fs;
sys.Ts = Ts;
sys.fmin = fmin;
sys.fmax = fmax;

for i = 1:input_dim
    u = u_time_domain{i};
    inject = inject_time_domain{i};
    for j = 1:output_dim
        y = y_time_domain{i, j};
        [h, coh, freqs] = bode_tfest_welch2(inject, u, y, sys.Fs, sys.fmin, sys.fmax);
        sys.tf{i, j}  = h;
        sys.coh{i, j} = coh;
    end
end
sys.freqs = freqs;

for i = 1:input_dim
    sys.siso_plant_frd_cont{i} = frd(sys.tf{i, i}, sys.freqs);
    sys.siso_plant_frd_disc{i} = frd(sys.tf{i, i}, sys.freqs, sys.Ts);
end

%% 绘制SISO传递函数
for i = 1:input_dim
    h = sys.tf{i, i};
    coh = sys.coh{i, i};
    freqs = sys.freqs;
    displayFreqdomain(h, coh, freqs, ...
        "Title", input_dim_name{i} + "->" + output_dim_name{i});
end

%% 绘制MIMO传递函数
displayComplexTransferMatrix(sys.tf, sys.freqs, input_dim_name, output_dim_name);

%% 绘制SISO惯性量
sys.mass = cell(input_dim,1);
sys.mass_mean = cell(input_dim,1);
sys.f1 = 5; sys.f2 = 14;

for i = 1:input_dim
    sys.mass{i} = 1 ./ ((2*pi*sys.freqs).^2 .*abs(sys.tf{i,i}));
end

for i = 1:input_dim
    sys.mass_mean{i} = mean(sys.mass{i}((sys.f1 <= sys.freqs) & (sys.freqs <= sys.f2)));
end


%% 函数定义
function displayTimedomain(u, y, varargin)
    % 1. 设置输入解析器
    p = inputParser;
    
    % 定义必填参数
    addRequired(p, 'u', @(x) true);
    addRequired(p, 'y', @(x) true);
    
    % 2. 定义带缺省值的可选参数 (使用名称-值对的方式)
    defaultFigid = 0;
    defaultTitle = '';
    defaultXLabel = "samples";
    defaultYLabel = '';
    defaultTs = -1;
    
    addParameter(p, 'Fig_id', defaultFigid);
    addParameter(p, 'Title1', defaultTitle);
    addParameter(p, 'Title2', defaultTitle);
    addParameter(p, 'YLabel1', defaultYLabel);
    addParameter(p, 'YLabel2', defaultYLabel);
    addParameter(p, 'XLabel', defaultXLabel);
    addParameter(p, 'Ts', defaultTs);
    
    % 3. 解析传入的参数
    parse(p, u, y, varargin{:});
    
    % 4. 提取解析后的结果
    u = p.Results.u;
    y = p.Results.y;

    fig_id = p.Results.Fig_id;
    titleStr1 = p.Results.Title1;
    titleStr2 = p.Results.Title2;
    yLabelStr1 = p.Results.YLabel1;
    yLabelStr2 = p.Results.YLabel2;
    xLabelStr = p.Results.XLabel;
    Ts = p.Results.Ts;

    if (Ts > 0)
        t1 = (0:length(u)-1) * Ts;
        t2 = (0:length(y)-1) * Ts;
    else
        t1 = (1:length(u));
        t2 = (1:length(y));
    end

    if ((xLabelStr == defaultXLabel) && Ts > 0)
        xLabelStr = "time/s";
    end
    
    % 5. 绘图逻辑
    if fig_id > 0
        figure(fig_id);
    else
        figure;
    end
    clf; % 清除当前图形窗口的旧内容
    
    subplot(2,1,1);
    plot(t1, u);
    xlabel(xLabelStr); % 应用横坐标标签
    ylabel(yLabelStr1); % 应用纵坐标标签
    title(titleStr1);   
    grid on;
    
    subplot(2,1,2);
    plot(t2, y);
    xlabel(xLabelStr); % 应用横坐标标签
    ylabel(yLabelStr2); 
    title(titleStr2);
    grid on;
end

function displayFreqdomain(h, coh, freqs, varargin)
    % 1. 设置输入解析器
    p = inputParser;
    
    % 定义必填参数
    addRequired(p, 'h', @(x) true);
    addRequired(p, 'coh', @(x) true);
    addRequired(p, 'freqs', @(x) true);

    % 2. 定义带缺省值的可选参数 (使用名称-值对的方式)
    defaultFigid = 0;
    defaultTitle = '';
    defaultXLabel = "frequency[Hz]";
    defaultYLabel1 = 'Magnitude[dB]';
    defaultYLabel2 = 'Phase[°]';
    defaultYLabel3 = 'Coherence';
    
    addParameter(p, 'Fig_id', defaultFigid);
    addParameter(p, 'Title', defaultTitle);
    addParameter(p, 'YLabel1', defaultYLabel1);
    addParameter(p, 'YLabel2', defaultYLabel2);
    addParameter(p, 'YLabel3', defaultYLabel3);
    addParameter(p, 'XLabel', defaultXLabel);

    % 3. 解析传入的参数
    parse(p, h, coh, freqs, varargin{:});
    
    % 4. 提取解析后的结果
    h = p.Results.h;
    coh = p.Results.coh;
    freqs = p.Results.freqs;

    fig_id = p.Results.Fig_id;
    titleStr = p.Results.Title;
    yLabelStr1 = p.Results.YLabel1;
    yLabelStr2 = p.Results.YLabel2;
    yLabelStr3 = p.Results.YLabel3;
    xLabelStr = p.Results.XLabel;

    if fig_id > 0
        figure(fig_id);
    else
        figure;
    end
    clf;
    subplot(3,1,1);
    semilogx(freqs, 20*log10(abs(h)));
    xlabel(xLabelStr); ylabel(yLabelStr1);
    title(titleStr);
    grid on;
    
    subplot(3,1,2);
    semilogx(freqs, rad2deg(unwrap(angle(h))));
    xlabel(xLabelStr); ylabel(yLabelStr2);
    grid on;
    
    subplot(3,1,3);
    semilogx(freqs, coh);
    xlabel(xLabelStr); ylabel(yLabelStr3);
    grid on;
end

function displayComplexTransferMatrix(H, freqs, input_dim_name, output_dim_name)
    [M, N, ~] = size(H);

    figure;
    for m = 1:M
        for n = 1:N
            i = (m-1)*N + n;
            subplot(M,N,i);
            yyaxis left;  semilogx(freqs, 20*log10(abs(squeeze(H{m,n})))); 
            yyaxis right; semilogx(freqs, rad2deg(unwrap(angle(squeeze(H{m,n})))));
            title_str = sprintf("%s->%s", input_dim_name{m}, output_dim_name{n});
            title(title_str); hold on;
            grid on;
        end
    end
end