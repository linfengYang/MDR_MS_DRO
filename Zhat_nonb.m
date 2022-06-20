%%  Zhat_nonb  by yy 2022.4.10
function [Z] = Zhat_nonb(Y,t,T,W,n,C,bN)

Z = [];
for tt = 1:t
    if tt == t
        Z = [Z,0,Y(n,t) - C];
    else
        Z = [Z,0,0];
    end
end

for tt = 1:t
    for w = 1:W
        Z = [Z,0,Y(n,T+1+(Sumt(t-1)+(tt-1))*W*bN+(w-1)*bN :T+(Sumt(t-1)+(tt-1))*W*bN+w*bN)];
    end
end
   
end