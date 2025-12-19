function [Numberindex] = findindex(rowindex,H_rowindex,colindex,mod)
if mod == 1
    len = 4;
elseif mod == 2
    len = 2;
else
    error('mod只能输入为1或2！');
end
for i = 1:len
    if colindex == H_rowindex(rowindex,i)
        Numberindex = i;
        break;
    end
end