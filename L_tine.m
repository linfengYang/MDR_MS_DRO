%%  Ltine  by yy 2021.8.29
function [ksi_tine] = L_tine(ksi,bN,bp)

ksi_tine = [];
if (bN == 1)
    ksi_tine = ksi;
else
    for i = 1:bN-1 %Ö¸Ê¾º¯Êý
        if (ksi < bp(i+1,1))
            ksi_tine = [ksi_tine; 0];
        else
            ksi_tine = [ksi_tine; 1];
        end
    end
end
    

end