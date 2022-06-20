%%  vers  by yy 2021.8.29
function [vers] = vers(bN,bp,bpr,bp0)

vers = [];
bp = [bp0; bp; bpr];
for j = 1:bN+1
    mid = [];
    mid = [mid; bp(j,1)];
    mid = [mid; L_hat(bp(j,1),bN,bp)];
    mid = [mid; L_tine(bp(j,1),bN,bp)];
    vers = [vers, mid];
end
%极限情况
for j = 2:bN
    mid = [];
    mid = [mid; bp(j,1)];
    mid = [mid; L_hat(bp(j,1),bN,bp)];
    mid = [mid; L_tine(bp(j,1)-0.001,bN,bp)];
    vers = [vers, mid];
end

end