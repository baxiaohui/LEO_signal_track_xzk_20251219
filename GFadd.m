function [GFadd_res] = GFadd(addend, add, m)
% 函数功能：实现伽罗华域的加法
% Input：
%   addend：被加数，向量法表示
%   add：加数，向量法表示
%   m：阶数
% Output：
%   GFadd_res：加法结果，向量法表示
len = length(addend);
ele = 2.^[m-1:-1:0];
addend_bin = deci2bina(addend, m);
add_bin = deci2bina(add, m);
for i = 1:len
    GFadd_res_bin = mod(addend_bin(i,:) + add_bin, 2);
    GFadd_res(i,1) = GFadd_res_bin * ele';
end