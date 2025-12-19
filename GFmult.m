function [GFmult_res] = GFmult(GFmultiplicand, GFmultiplier, v2p_tab, p2v_tab, m)
% 函数功能：完成一个域元素矩阵或一个域元素向量与一个域元素的乘法
% Input：
%   GFmultiplicand：被乘数，是一个域元素矩阵或域元素向量
%   GFmultiplier：乘数，是一个域元素
%   v2p_tab：向量表示法到幂次表示法的映射
%   p2v_tab：幂次表示法到向量表示法的映射
%   m：幂次的阶数
% Output：
%   GFmult_res：乘法结果,向量法表示
GFq = 2^m - 1;
ele = 2.^[m-1:-1:0];
[M,N] = size(GFmultiplicand);
GFmult_res = zeros(M,N);
GFmultiplier_pow = v2p_tab(GFmultiplier + 1);
if GFmultiplier_pow == -1
    GFmult_res = zeros(M,N);
else
    for i = 1:M
        for j = 1:N
            GFmultiplicand_pow = v2p_tab(GFmultiplicand(i,j)+1);
            if GFmultiplicand_pow == -1
                GFmult_res(i,j) = 0;
            else
                tmp = mod(GFmultiplicand_pow + GFmultiplier_pow, GFq);
                GFmult_res(i,j) = p2v_tab(tmp + 2,:) * ele';
            end
        end
    end
end