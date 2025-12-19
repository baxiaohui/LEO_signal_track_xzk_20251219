function GF_domain = GF_gene(m)
% 本函数用于生成 2^m 的伽罗华域表，其中 3≤m≤10
% GF_domain是输出的伽罗华域表，共 2^m 组，每组长为 m 位

if m<3 || m>10
    error('m只能在3到10之间！');
else
    X_M = [1 1 0 0 0 0 0 0 0 0;...
           1 1 0 0 0 0 0 0 0 0;...
           1 0 1 0 0 0 0 0 0 0;...
           1 1 0 0 0 0 0 0 0 0;...
           1 0 0 1 0 0 0 0 0 0;...
           1 0 1 1 1 0 0 0 0 0;...
           1 0 0 0 1 0 0 0 0 0;...
           1 0 0 1 0 0 0 0 0 0;];
end
x_m = X_M(m-2, 1:m);
block_num = 2^m;
GF_domain = zeros(block_num, m);
GF_domain(2, :) = eye(1, m);
flag = 0;

for i = 3 : block_num
    row_tmp = GF_domain(i-1, :);
    flag = row_tmp(end);
    if flag == 0
        GF_domain(i, :) = [row_tmp(end),row_tmp(1:end-1)];
    elseif flag == 1
        GF_domain_tmp = [0, row_tmp(1:end-1)] + x_m;
        GF_domain(i, :) = mod(GF_domain_tmp, 2);
    end
end