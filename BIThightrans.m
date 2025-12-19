function dec_values = BIThightrans(bin_array)
    % 确保输入是行向量
    bin_array = bin_array(:)';
    
    % 检查输入长度是否为6的倍数
    if mod(length(bin_array), 6) ~= 0
        error('输入二进制数组的长度必须是6的倍数。');
    end
    
    % 将二进制数组重塑为6行N列的矩阵
    groups = reshape(bin_array, 6, []);
    
    % 定义权值向量（2^5 到 2^0）
    weights = 2.^(5:-1:0);
    
    % 计算每组的十进制值
    dec_values = weights * groups;
    
    % 转换为行向量输出
    dec_values = dec_values(:)';
end