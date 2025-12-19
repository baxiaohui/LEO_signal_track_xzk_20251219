function B_Str = deci2bina(Z,n)
% This function is used for converting a nonnegtive integer into a
% n-bit vector that can represent this integer
Z_len = length(Z);
B_Str = zeros(Z_len,n);
for i = 1:Z_len
    for j = 1:n
        if (Z(i) >= 2^(n-j))
            B_Str(i,j) = 1;
            Z(i) = Z(i) - 2^(n-j);
        end
    end
end