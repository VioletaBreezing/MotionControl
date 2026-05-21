function data = readSingleColumnData(filename)
    % 读取单列浮点数TXT文件的函数
    % 输入参数:
    %   filename - 文件名（字符串），包含路径信息
    % 输出参数:
    %   data - 包含文件中所有浮点数的列向量
    
    % 方法1: 使用load函数（适用于纯数值文件）
    try
        data = load(filename);
        % 如果文件只有一列，load返回的是列向量；否则转换为列向量
        if size(data, 2) > 1
            data = data(:);
        end
    catch
        % 如果load失败，则使用textscan方法
        fileID = fopen(filename, 'r');
        if fileID == -1
            error('无法打开文件: %s', filename);
        end
        
        % 读取每一行的浮点数
        data = textscan(fileID, '%f', 'Delimiter', '\n');
        data = data{1};  % 提取cell数组中的数值向量
        fclose(fileID);
    end
    
    fprintf('成功读取 %d 个数据点\n', length(data));
end
