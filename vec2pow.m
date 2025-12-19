function [v2m] = vec2pow(p2v, m)
% 函数功能：实现域元素的向量表示法到幂表示法的查找表
% Input：
%   p2v：幂表示法到向量表示法的查找表
%   m：幂的阶数
% Output：
%   v2m：向量表示法到幂表示法的查找表

ele = 2.^[m-1:-1:0];
bina2deci = p2v * ele';
[~,v2m_tmp] = sort(bina2deci,1,'ascend');
v2m = v2m_tmp - 2;