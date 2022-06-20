%%  Zbin  by yy 2022.3.28
function [Z] = Zbin(Y,t,T,W,n,C,bpN,bN)

Z = [];
for tt = 1:t
    if tt==t
        Z = [Z,0,0,Y(n,t) - C];
    else
        Z = [Z,0,0,0];
    end
end
for tt = 1:t
    for w = 1:W
        Z = [Z,0,zeros(1,bN),Y(n,T+1+(Sumt(t-1)+(tt-1))*W*bpN+(w-1)*bpN :T+(Sumt(t-1)+(tt-1))*W*bpN+w*bpN)];
    end
end
   
end

% function [Z] = Zbin(Y,t,T,W,n,C,bpN,bN)
% 
% Z = [0,0,Y(n,t) - C];
% for tt = 1:t
%     for w = 1:W
%         Z = [Z,0,zeros(1,bN),Y(n,T+1+(Sumt(t-1)+(tt-1))*W*bpN+(w-1)*bpN :T+(Sumt(t-1)+(tt-1))*W*bpN+w*bpN)];
%     end
% end
%    
% end