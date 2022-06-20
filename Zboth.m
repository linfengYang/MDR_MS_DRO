%%  Zboth  by yy 2022.3.25
function [Z] = Zboth(Y1,Y2,t,T,W,n,C1,C2,bpN,bN)

Z = [0,Y1(n,t) - C1,Y2(n,t) - C2];
for tt = 1:t
    for w = 1:W
        Z = [Z,0,Y1(n,T+1+(Sumt(t-1)+(tt-1))*W*bN+(w-1)*bN :T+(Sumt(t-1)+(tt-1))*W*bN+w*bN),zeros(1,bpN)];
    end
end
    

end