%%  U_matrix  by yy 2021.10.29
% 考虑写一般形式循环太多，这里写为分段数为3的情况
function [U] = U_matrix(verts_noksi,t)

U = [];
v = [];
[row,col] = size(verts_noksi);
I = row*2;
for i = 1:I % i的个数等于风电机组数*（分段数加结点数）
    if i <= row
        v(i,:) = verts_noksi(i,1+(t-1)*8:4+(t-1)*8);
    else
        v(i,:) = verts_noksi(i-row,5+(t-1)*8:8*t);
    end
end
for i = 1:I
    U = blkdiag(U,v(i,:));
end
% U = blkdiag(v(1,:),v(2,:),v(3,:),v(4,:),v(5,:),v(6,:),v(7,:),v(8,:),v(9,:),v(10,:));
% 如果写成一般形式，这里可以用U = blkdiag(U,v(i,:))循环实现

end