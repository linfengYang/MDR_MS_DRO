%%  Lhat  by yy 2021.8.27
% function [ksi_hat] = L_hat(ksi,bN,bp)
% 
% ksi_hat = [];
% if (bN == 1)
%     ksi_hat = ksi;
% else
%     for i = 1:bN
%         if (i == 1)
%             mid = min(ksi,bp(i+1,1));
%             ksi_hat = [ksi_hat;mid];
%         else
%             mid = min(ksi,bp(i+1,1));
%             mid = mid - bp(i,1);
%             mid = max(0,mid);
%             ksi_hat = [ksi_hat;mid];
%         end       
%     end
% end
% 
% end
function [ksi_hat] = L_hat(ksi,bN,bp)

ksi_hat = [];
if (bN == 1)
    ksi_hat = ksi;
else
    for i = 1:bN
        mid = min(ksi,bp(i+1,1));
        mid = mid - bp(i,1);
        mid = max(0,mid);
        ksi_hat = [ksi_hat;mid];
    end
end

end