%%  Zhat  by yy 2022.3.25
function [Z] = Zhat(Y,t,T,W,n,C,bpN,bN)

Z = [];
for tt = 1:t
    if tt == t
        Z = [Z,0,Y(n,t) - C,0];
    else
        Z = [Z,0,0,0];
    end
end

for tt = 1:t
    for w = 1:W
        Z = [Z,0,Y(n,T+1+(Sumt(t-1)+(tt-1))*W*bN+(w-1)*bN :T+(Sumt(t-1)+(tt-1))*W*bN+w*bN),zeros(1,bpN)];
    end
end
   
end

% function [Z] = Zhat(Y,t,T,W,n,C,bpN,bN)
% if bpN == 0
%     Z = [0,Y(n,t) - C];
% else
%     Z = [0,Y(n,t) - C,0];
% end
% for tt = 1:t
%     for w = 1:W
%         Z = [Z,0,Y(n,T+1+(Sumt(t-1)+(tt-1))*W*bN+(w-1)*bN :T+(Sumt(t-1)+(tt-1))*W*bN+w*bN),zeros(1,bpN)];
%     end
% end
%    
% end
