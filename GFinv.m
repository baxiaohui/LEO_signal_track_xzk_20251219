function [GFinv_res] = GFinv(GFele, v2p_tab, p2v_tab, m)
% 函数功能：求出一个域元素的逆
% Input：
%   GFele：待求逆的域元素
%   v2p_tab：向量表示法到幂次表示法的查找表
%   p2v_tab：幂次表示法到向量表示法的查找表
%   m：幂次的阶数
% Output：
%   GFinv_res：求逆的结果
ele = 2.^[m-1:-1:0];
GFnum = 2^m - 1;
GFele_pow = v2p_tab(GFele+1);
GFele_pow_inv = mod(GFnum - GFele_pow, GFnum);
GFele_vec_inv = p2v_tab(GFele_pow_inv+2,:);
GFinv_res = GFele_vec_inv * ele';