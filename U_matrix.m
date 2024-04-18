%%  U_matrix  by yy 2021.10.29
function [U] = U_matrix(verts_noksi,t)

U = [];
v = [];
[row,col] = size(verts_noksi);
I = row*2;
for i = 1:I 
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
end
