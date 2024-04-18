%%  BDR_UCC_DRO_Last  by yy 2022.04.10

FileName_SCUC = 'data/SCUC6.txt';
FileName_Wind = 'data\Wind_power.xlsx';

SCUC_data = ReadDataSCUC(FileName_SCUC);
Wind_data = ReadWindData(SCUC_data,FileName_Wind);

% T = SCUC_data.totalLoad.T;  % period T
T = 12;  
G = SCUC_data.units.N;      % number of thermal units
N = SCUC_data.baseparameters.busN;  % number of buses
W = SCUC_data.Wind.wind_Number; % number of wind units
M = size(Wind_data{1,1},2); % number of samples
Gbus = SCUC_data.units.bus_G; % thermal units
Wbus = SCUC_data.Wind.wind_Node; % wind units
Dbus = SCUC_data.busLoad.bus_PDQR; % load buses
Pramp = SCUC_data.units.ramp; %ramp constraints
Pstart = 0.5 * SCUC_data.units.PG_up; %startup ramp constraints
Pshut = SCUC_data.units.PG_up; %shutdown ramp constraints 
Pup = SCUC_data.units.PG_up; %upper bound of generation
Plow = SCUC_data.units.PG_low; %lower bound
Ton = SCUC_data.units.T_on; %min on time
Toff = SCUC_data.units.T_off; %min off time
L = 4; %parameter for linearize
fa = SCUC_data.units.alpha; %constant coefficient for generation cost
fb = SCUC_data.units.beta; %linear term coefficient for generation cost
fc = SCUC_data.units.gamma; % quadratic term for generation cost
Ccold = SCUC_data.units.start_cost*2; % cost of cold start
Chot = SCUC_data.units.start_cost; % cost of hot start
Tcold = min(max(Toff,1),T); %time of cold start
Lup = 0.3*ones(N,1); % upper bound of load loss
Llow = 0*ones(N,1); %lower bound of load loss

% construct coefficient matrix of DC network B
type_of_pf = 'DC';
Y = SCUC_nodeY(SCUC_data,type_of_pf);
B = -Y.B;

u0 = zeros(G,1); %initial state of units，assume all shutdown at beggin
T0 = zeros(G,1); %min periods that unit should keep in initial state
P0 = zeros(G,1); %intial generation
U0 = zeros(G,1); 
L0 = zeros(G,1); 
tao0 = Toff+Tcold+1; 

for index = 1:G
    mid_u = min(T,u0(index) * (Ton(index) - T0(index)));
    mid_l = min(T,(1-u0(index)) * (Toff(index) + T0(index)));
    U0(index) = max(0,mid_u);
    L0(index) = max(0,mid_l);
end

all_branch.I = [ SCUC_data.branch.I; SCUC_data.branchTransformer.I ]; % all branch begin bus
all_branch.J = [ SCUC_data.branch.J; SCUC_data.branchTransformer.J ]; % all branch end bus
all_branch.P = [ SCUC_data.branch.P; SCUC_data.branchTransformer.P ]; % transmission limitation of each branch

% calculate uncertainty set
PW = []; % put wind generation variables in a matrix
for i = 1:W
    PW = [PW; Wind_data{1,i}];
end
PW = [PW(1:T,:);PW(24+1:24+T,:)];
% for i = 1:W
%     PW = [PW; 0.05*ones(24,30)];
% end
% average value of wind generation
PW_miu = zeros(W*T,1);
for i = 1:M
    PW_miu = PW_miu + PW(:,i);
end
PW_miu = PW_miu / M;
% error of wind generation
VW = zeros(W*T,M);
for j = 1:M
    VW(:,j) =  PW(:,j) - PW_miu;
end
VW_positive =  max(VW,[],2);
VW_negative =  min(VW,[],2);
% ksi support set
ksi_positive = [1;VW_positive];
ksi_negative = [1;VW_negative];
% end of calculate uncertainty set

% non-linear decision rule
bN = 3; % number of segments
bpN = bN-1; % inflection points
bp = zeros(bN - 1,1); % equally spaced vector
kl = 1+bpN+bN; % dimension of one ksi
verts = []; % set of vertices
for t = 1:T
    for i = 1:W
        step = (VW_positive((t-1)*W+i) - VW_negative((t-1)*W+i)) / bN; %calculate the length of each segment
        bp = linspace(VW_negative((t-1)*W+i)+ step,VW_positive((t-1)*W+i)-step,bN-1)';
        verts = [verts, vers(bN,bp,VW_positive((t-1)*W+i),VW_negative((t-1)*W+i))];
    end
end
verts_nonksi = verts(2:6,:);
vN = length(verts)/(W*T); %number of verts
% verts = [ones(2*bN,(bN+1)*W), verts];
% end of non-linear decision rule

% use decision rules handle all ksi vector
for m = 1:M
    ksi = VW(:,m);
    ksi_hat = [];
    ksi_tine = [];
    ksi_hat2 = [];
    ksi_tine2 = [];
    ksi_wan_now = [];
    for t = 1:T
        for i = 1:W
            step = (VW_positive((t-1)*W+i) - VW_negative((t-1)*W+i)) / bN;
            bp = linspace(VW_negative((t-1)*W+i)+ step,VW_positive((t-1)*W+i)-step,bN-1)';
            bp = [VW_negative((t-1)*W+i);bp;VW_positive((t-1)*W+i)];
            ksi_hat = L_hat(ksi((t-1)*W+i),bN,bp);
            ksi_tine = L_tine(ksi((t-1)*W+i),bN,bp);
            ksi_hat2 = [ksi_hat2; L_hat(ksi((i-1)*T+t),bN,bp)];
            ksi_tine2 = [ksi_tine2;  L_tine(ksi((i-1)*T+t),bN,bp)];
            ksi_wan_now = [ksi_wan_now;ksi((t-1)*W+i);ksi_hat;ksi_tine];
        end
    end
%     ksi_wan{m} = [ones(T,1); ones(T,1); ksi_hat; ksi_tine]; %arrange as this order :ksi_hat0;ksi_tine0;ksi_hat;ksi_tine
    ksi_wan{m} = [ ones(T,1);ones(T,1);ones(T,1);ksi_wan_now]; %arrange as this order: ksi0;ksi_hat0;ksi_tine0;ksi;ksi_hat;ksi_tine
    ksi_wan2{m} = [ ksi_hat2; ksi_tine2]; %arrange as: ksi_hat;ksi_tine
end
% end of use decision rules handle all ksi vector

% % construct wasserstein ambiguity set
% eta = 0.95; % confident level
% rho = sdpvar(1); %radius of ambiguity set
% sum_C = 0;
% for i = 1:M
%     mid = PW(:,i) - PW_miu;
%     mid = rho * norm(mid,1)^2;
%     mid = exp(mid);
%     sum_C = sum_C + mid;
% end
% sum_C = sum_C / M;
% obj_C = 1 + log(sum_C);
% Constraints_C =  rho >= 0;
% Objective_C = 2 *  ( ((1 / (2 * rho)) * obj_C) ^(1/2) );
% options_C = sdpsettings('verbose',0,'debug',1,'savesolveroutput',1);%, 'fmincon.TolX',1e-4
% sol_C = optimize(Constraints_C,Objective_C,options_C);
% C = sol_C.solveroutput.fmin;
% mid_eps = log( 1 / (1-eta));
% mid_eps = mid_eps / M;
% eps = C * sqrt(mid_eps);
% % end of construct wasserstein ambiguity set
eps = 4.5;
% construct closed convex hull
% verts already finished in last section

% end of construct closed convex hull

%calculate ksi length
leng_k = 0;
for k = 1:T
    leng_k = leng_k+k;
end
%end of calculate ksi length

% construct continuous constraints
% define continuous variables 
YG = sdpvar(G,T+W*bN*leng_k); %generation decision variable for thermal units
% YG = 1:64;
% YG = [YG;YG;YG;YG;YG;YG];
Ytheta = sdpvar(N,T+W*bN*leng_k); %nodal phase angle decision variables
Yz = sdpvar(G,T+W*bN*leng_k); % cost variable for thermal units
YS = sdpvar(G,T+W*bN*leng_k); % startup cost variable for thermal units
YPw = sdpvar(W,T+W*bN*leng_k); % wind generation decision variable
Yl = sdpvar(N,T+W*bN*leng_k); %load loss decision variable
% integer variable
% not binary variable
Xu = intvar(G,T+W*bpN*leng_k); %unit state varibale
Xs = intvar(G,T+W*bpN*leng_k); %the action of turn on
Xd = intvar(G,T+W*bpN*leng_k); %the action of turn off

% construct constraints
cons = [];
cons = [cons, -1 <= Xu <= 1]; % integer variable's range
cons = [cons, -1 <= Xs <= 1];
cons = [cons, -1 <= Xd <= 1];

% reference bus
Ytheta(1,:) = zeros(1,T+W*bN*leng_k);
% end of reference bus

% keep in initial state
for g = 1:G 
    for t = 1:(U0(g)+L0(g))
        Xu(g,T+1+Sumt(t-1)*W*bpN:T+Sumt(t)*W*bpN) = zeros(1,t*W*bpN);
        Xu(g,t) = 0;
    end
end
% end of keep in initial state

% construct inequality constraints
Uc_l = []; %matrix of verts
Uc_c = []; %lamda coefficient matrix
Uc = [];
Uconstant_0 = []; 
Uconstant_0 = blkdiag(1, Uconstant_0);
Uconstant_0 = blkdiag(1, Uconstant_0);
Uconstant_0 = blkdiag(1, Uconstant_0);
Ulamda_0 = [];
Ulamda_0 = blkdiag(-1, Ulamda_0); 
Ulamda_0 = blkdiag(-1, Ulamda_0); 
Ulamda_0 = blkdiag(-1, Ulamda_0); 
for t = 1:T
    Uconstant = []; 
    for w = 1:W
        Uconstant = blkdiag(Uconstant,ones(1,vN));
    end
    Ulamda = blkdiag(-verts(:,1+(t-1)*2*vN:vN+(t-1)*2*vN),-verts(:,1+vN+(t-1)*2*vN:2*vN*t)); 
    Uc_l = blkdiag(Ulamda_0,Uc_l); 
    Uc_l = blkdiag(Uc_l,Ulamda); 
    Uc_c = blkdiag(Uconstant_0,Uc_c); 
    Uc_c = blkdiag(Uc_c,Uconstant); 
    Uc = [Uc_l;Uc_c]; 
    Wc_l = diag(ones(t*W*kl+3*t,1)); 
    Wc_c = zeros(t*W+3*t,t*W*kl+3*t);
    Wc = [Wc_l;Wc_c];
    hc = [zeros(t*W*kl + 3*t,1);ones(t*W+3*t,1)]; 
    for n = 1:N % traverse all buses
        % constraints for upper bound of load loss 
        Lamda_lup = [sdpvar(t*W*kl + 3*t,1); sdpvar(t*W + 3*t,1)];
        Lamda2_lup = sdpvar((t*W)*vN + 3*t,1); 
        Zlup = -Zhat(Yl,t,T,W,n,Lup(n,1),bpN,bN);
        cons = [cons, Wc'* Lamda_lup == Zlup'];
        cons = [cons, Uc'* Lamda_lup + Lamda2_lup == 0];
        cons = [cons, hc'* Lamda_lup >= 0];
        cons = [cons, Lamda2_lup >= 0];
        % constraints for lower bound of load loss 
        Lamda_llow = [sdpvar(t*W*kl + 3*t,1); sdpvar(t*W + 3*t,1)]; 
        Lamda2_llow = sdpvar((t*W)*vN + 3*t,1); 
        Zllow = Zhat(Yl,t,T,W,n,Llow(n,1),bpN,bN);
        cons = [cons, Wc'* Lamda_llow == Zllow'];
        cons = [cons, Uc'* Lamda_llow + Lamda2_llow == 0];
        cons = [cons, hc'* Lamda_llow >= 0];
        cons = [cons, Lamda2_llow >= 0];
        %check if current bus contains wind generation
        if ismember(n,Wbus)
            index = find(Wbus==n);
            Lamda_w = [sdpvar(t*W*kl + 3*t,1); sdpvar(t*W + 3*t,1)];
            Lamda2_w = sdpvar((t*W)*vN + 3*t,1);
            Aw = [PW_miu((t-1)*W+index,1)-YPw(index,t), -YPw(index,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN)];
            Zw = [];
            for tt = 1:t
                if tt == t
                    Zw = [Zw,0,Aw(1,1),0];
                else
                    Zw = [Zw,0,0,0];
                end
            end
            for tt = 1:t
                for w = 1:W
                    if (tt==t && w==index)
                        Zw = [Zw,1,Aw(1,2+(tt-1)*W*bN+(w-1)*bN:1+(tt-1)*W*bN+w*bN),zeros(1,bpN)];
                    else
                        Zw = [Zw,0,Aw(1,2+(tt-1)*W*bN+(w-1)*bN:1+(tt-1)*W*bN+w*bN),zeros(1,bpN)];
                    end
                end
            end
            cons = [cons, Wc'* Lamda_w == Zw'];
            cons = [cons, Uc'* Lamda_w + Lamda2_w == 0];
            cons = [cons, hc'* Lamda_w >= 0];
            cons = [cons, Lamda2_w >= 0];
            % lower bound of wind generation
            Lamda_wlow = [sdpvar(t*W*(1+bN+bpN) + 3*t,1); sdpvar(t*W + 3*t,1)]; 
            Lamda2_wlow = sdpvar((t*W)*vN + 3*t,1); 
            Zwlow = Zhat(YPw,t,T,W,index,0,bpN,bN);
            cons = [cons, Wc'* Lamda_wlow == Zwlow'];
            cons = [cons, Uc'* Lamda_wlow + Lamda2_wlow == 0];
            cons = [cons, hc'* Lamda_wlow >= 0];
            cons = [cons, Lamda2_wlow >= 0];
        end
        % check if current bus contains thermal unit
        if ismember(n,Gbus) 
            index = find(Gbus==n);
            %lowest generation cost 
            Lamda_wan = [sdpvar(t*W*kl + 3*t,1); sdpvar(t*W + 3*t,1)];
            Lamda2_wan = sdpvar((t*W)*vN + 3*t,1);
            Zwan = Zhat(YS,t,T,W,index,0,bpN,bN);
            cons = [cons, Wc'* Lamda_wan == Zwan'];
            cons = [cons, Uc'* Lamda_wan + Lamda2_wan == 0];
            cons = [cons, hc'* Lamda_wan >= 0];
            cons = [cons, Lamda2_wan >= 0];
            
            % linearize the cost of generation
            for l = 0:(L-1)
                Lamda_cost = [sdpvar(t*W*kl + 3*t,1); sdpvar(t*W + 3*t,1)];
                Lamda2_cost = sdpvar((t*W)*vN + 3*t,1);
                p_i_l = Plow(index) + (Pup(index) - Plow(index)) / L * l;
                Acost = [Yz(index,t), Yz(index,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN)] - (2*fc(index)*p_i_l+fb(index)) * [YG(index,t),YG(index,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN)];
                Bcost = -(fa(index)-fc(index)*p_i_l*p_i_l) * [Xu(index,t),Xu(index,T+1+Sumt(t-1)*W*bpN:T+Sumt(t)*W*bpN)];
                Zcost = [];
                for tt = 1:t
                    if tt == t
                        Zcost = [Zcost,0,Acost(1,1),Bcost(1,1)];
                    else
                        Zcost = [Zcost,0,0,0];
                    end
                end
                for tt = 1:t 
                    for w = 1:W
                        Zcost = [Zcost,0,Acost(1,2+(tt-1)*W*bN+(w-1)*bN:1+(tt-1)*W*bN+w*bN),Bcost(1,2+(tt-1)*W*bpN+(w-1)*bpN:1+(tt-1)*W*bpN+w*bpN)];
                    end
                end
                cons = [cons, Wc'* Lamda_cost == Zcost'];
                cons = [cons, Uc'* Lamda_cost + Lamda2_cost == 0];
                cons = [cons, hc'* Lamda_cost >= 0];
                cons = [cons, Lamda2_cost >= 0];
            end
            
            % cost of startup
            if (t-tao0(index)+1<=0)&&(max(0,-T0(index))<abs(t-tao0(index)+1))
                fit = 1;
            else
                fit = 0;
            end
            tao_it = max(1,t-Toff(index)-Tcold(index));
            Lamda_open = [sdpvar(t*W*kl + 3*t,1); sdpvar(t*W + 3*t,1)];
            Lamda2_open = sdpvar((t*W)*vN + 3*t,1);
            Ccost = Chot(index) - Ccold(index); 
            Aopen = [YS(index,t)-fit*Ccost,YS(index,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN)];
            Bopen = zeros(1,1+t*W*bpN);
            for y = tao_it : t-1 
                Bopen = Bopen - [Xd(index,y),Xd(index,T+1+Sumt(y-1)*W*bpN:T+Sumt(y)*W*bpN), zeros(1,(t-y)*W*bpN)];
            end
            Bopen = Bopen + [Xs(index,t),Xs(index,T+1+Sumt(t-1)*W*bpN:T+Sumt(t)*W*bpN)]; 
            Bopen = Ccost*Bopen;
            Zopen = [];
            for tt = 1:t
                if tt == t
                    Zopen = [Zopen,0,Aopen(1,1),Bopen(1,1)];
                else
                    Zopen = [Zopen,0,0,0];
                end
            end
            for tt = 1:t
                for w = 1:W
                    Zopen = [Zopen,0,Aopen(1,2+(tt-1)*W*bN+(w-1)*bN:1+(tt-1)*W*bN+w*bN),Bopen(1,2+(tt-1)*W*bpN+(w-1)*bpN:1+(tt-1)*W*bpN+w*bpN)];
                end
            end
            cons = [cons, Wc'* Lamda_open == Zopen'];
            cons = [cons, Uc'* Lamda_open + Lamda2_open == 0];
            cons = [cons, hc'* Lamda_open >= 0];
            cons = [cons, Lamda2_open >= 0]; 
            
            %upper bound of generation
            Lamda_up = [sdpvar(t*W*kl + 3*t,1); sdpvar(t*W + 3*t,1)];
            Lamda2_up = sdpvar((t*W)*vN + 3*t,1);
            Aup = [YG(index,t),YG(index,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN)]; 
            Bup = Pup(index) * [Xu(index,t),Xu(index,T+1+Sumt(t-1)*W*bpN:T+Sumt(t)*W*bpN)]; 
            Zup = [];
            for tt = 1:t
                if tt == t
                    Zup = [Zup,0,-Aup(1,1),Bup(1,1)];
                else
                    Zup = [Zup,0,0,0];
                end
            end
            for tt = 1:t 
                for w = 1:W
                    Zup = [Zup,0,-Aup(1,2+(tt-1)*W*bN+(w-1)*bN:1+(tt-1)*W*bN+w*bN),Bup(1,2+(tt-1)*W*bpN+(w-1)*bpN:1+(tt-1)*W*bpN+w*bpN)];
                end
            end
            cons = [cons, Wc'* Lamda_up == Zup'];
            cons = [cons, Uc'* Lamda_up + Lamda2_up == 0];
            cons = [cons, hc'* Lamda_up >= 0];
            cons = [cons, Lamda2_up >= 0]; 
            
            %lower bound of generation
            Lamda_low = [sdpvar(t*W*kl + 3*t,1); sdpvar(t*W + 3*t,1)];
            Lamda2_low = sdpvar((t*W)*vN + 3*t,1);
            Alow = [YG(index,t),YG(index,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN)];
            Blow = Plow(index) * [Xu(index,t),Xu(index,T+1+Sumt(t-1)*W*bpN:T+Sumt(t)*W*bpN)];
            Zlow = [];
            for tt = 1:t
                if tt == t
                    Zlow = [Zlow,0,Alow(1,1),-Blow(1,1)];
                else
                    Zlow = [Zlow,0,0,0];
                end
            end
            for tt = 1:t 
                for w = 1:W
                    Zlow = [Zlow,0,Alow(1,2+(tt-1)*W*bN+(w-1)*bN:1+(tt-1)*W*bN+w*bN),-Blow(1,2+(tt-1)*W*bpN+(w-1)*bpN:1+(tt-1)*W*bpN+w*bpN)];
                end
            end
            cons = [cons, Wc'* Lamda_low == Zlow'];
            cons = [cons, Uc'* Lamda_low + Lamda2_low == 0];
            cons = [cons, hc'* Lamda_low >= 0];
            cons = [cons, Lamda2_low >= 0];
            
            %up ramping
            Lamda_ramp_up = [sdpvar(t*W*kl + 3*t,1); sdpvar(t*W + 3*t,1)];
            Lamda2_ramp_up = sdpvar((t*W)*vN + 3*t,1);
            if t == 1
                Aramp_up =  -[YG(index,t),YG(index,T+1:T+W*bN)]; %当t=1时，只剩s和p1两项
                Bramp_up = Pstart(index) * [Xs(index,t),Xs(index,T+1:T+W*bpN)];
            else
                Aramp_up = [YG(index,t-1),YG(index,T+1+Sumt(t-2)*W*bN:T+Sumt(t-1)*W*bN),zeros(1,W*bN)] - [YG(index,t),YG(index,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN)]; 
                Bramp_up = Pramp(index) * [Xu(index,t-1),Xu(index,T+1+Sumt(t-2)*W*bpN:T+Sumt(t-1)*W*bpN),zeros(1,W*bpN)] + Pstart(index) * [Xs(index,t),Xs(index,T+1+Sumt(t-1)*W*bpN:T+Sumt(t)*W*bpN)]; 
            end
            Zramp_up = [];
            for tt = 1:t
                if tt == t
                    Zramp_up = [Zramp_up,0,Aramp_up(1,1),Bramp_up(1,1)];
                else
                    Zramp_up = [Zramp_up,0,0,0];
                end
            end
            for tt = 1:t 
                for w = 1:W
                    Zramp_up = [Zramp_up,0,Aramp_up(1,2+(tt-1)*W*bN+(w-1)*bN:1+(tt-1)*W*bN+w*bN),Bramp_up(1,2+(tt-1)*W*bpN+(w-1)*bpN:1+(tt-1)*W*bpN+w*bpN)];
                end
            end
            cons = [cons, Wc'* Lamda_ramp_up == Zramp_up'];
            cons = [cons, Uc'* Lamda_ramp_up + Lamda2_ramp_up == 0];
            cons = [cons, hc'* Lamda_ramp_up >= 0];
            cons = [cons, Lamda2_ramp_up >= 0];
            
            %down ramping
            Lamda_ramp_down = [sdpvar(t*W*kl + 3*t,1); sdpvar(t*W + 3*t,1)];
            Lamda2_ramp_down = sdpvar((t*W)*vN + 3*t,1);
            if t == 1
                Aramp_down = [YG(index,t),YG(index,T+1:T+W*bN)]; 
                Bramp_down = Pshut(index) * [Xd(index,t),Xd(index,T+1:T+W*bpN)] + Pramp(index) * [Xu(index,t),Xu(index,T+1:T+W*bpN)];
            else
                Aramp_down = [YG(index,t),YG(index,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN)] - [YG(index,t-1),YG(index,T+1+Sumt(t-2)*W*bN:T+Sumt(t-1)*W*bN), zeros(1,W*bN)]; 
                Bramp_down = Pshut(index) * [Xd(index,t),Xd(index,T+1+Sumt(t-1)*W*bpN:T+Sumt(t)*W*bpN)] + Pramp(index) * [Xu(index,t),Xu(index,T+1+Sumt(t-1)*W*bpN:T+Sumt(t)*W*bpN)]; 
            end
            Zramp_down = [];
            for tt = 1:t
                if tt == t
                    Zramp_down = [Zramp_down,0,Aramp_down(1,1),Bramp_down(1,1)];
                else
                    Zramp_down = [Zramp_down,0,0,0];
                end
            end
            for tt = 1:t 
                for w = 1:W
                    Zramp_down = [Zramp_down,0,Aramp_down(1,2+(tt-1)*W*bN+(w-1)*bN:1+(tt-1)*W*bN+w*bN),Bramp_down(1,2+(tt-1)*W*bpN+(w-1)*bpN:1+(tt-1)*W*bpN+w*bpN)];
                end
            end
            cons = [cons, Wc'* Lamda_ramp_down == Zramp_down'];
            cons = [cons, Uc'* Lamda_ramp_down + Lamda2_ramp_down== 0];
            cons = [cons, hc'* Lamda_ramp_down >= 0];
            cons = [cons, Lamda2_ramp_down >= 0];
            
            %only binary
            %startpu u>0 constraints
            Lamda_u0 = [sdpvar(t*W*kl + 3*t,1); sdpvar(t*W + 3*t,1)];
            Lamda2_u0 = sdpvar((t*W)*vN + 3*t,1);
            Zu0 = Zbin(Xu,t,T,W,index,0,bpN,bN);
            cons = [cons, Wc'* Lamda_u0 == Zu0'];
            cons = [cons, Uc'* Lamda_u0 + Lamda2_u0== 0];
            cons = [cons, hc'* Lamda_u0 >= 0];
            cons = [cons, Lamda2_u0 >= 0];
            %u<1
            Lamda_u1 = [sdpvar(t*W*kl + 3*t,1); sdpvar(t*W + 3*t,1)];
            Lamda2_u1 = sdpvar((t*W)*vN + 3*t,1);
            Zu1 = -Zbin(Xu,t,T,W,index,1,bpN,bN);
            cons = [cons, Wc'* Lamda_u1 == Zu1'];
            cons = [cons, Uc'* Lamda_u1 + Lamda2_u1== 0];
            cons = [cons, hc'* Lamda_u1 >= 0];
            cons = [cons, Lamda2_u1 >= 0];
            %s>0 constaints
            Lamda_s0 = [sdpvar(t*W*kl + 3*t,1); sdpvar(t*W + 3*t,1)];
            Lamda2_s0 = sdpvar((t*W)*vN + 3*t,1);
            Zs0 = Zbin(Xs,t,T,W,index,0,bpN,bN);
            cons = [cons, Wc'* Lamda_s0 == Zs0'];
            cons = [cons, Uc'* Lamda_s0 + Lamda2_s0== 0];
            cons = [cons, hc'* Lamda_s0 >= 0];
            cons = [cons, Lamda2_s0 >= 0];
            %s<1
            Lamda_s1 = [sdpvar(t*W*kl + 3*t,1); sdpvar(t*W + 3*t,1)];
            Lamda2_s1 = sdpvar((t*W)*vN + 3*t,1);
            Zs1 = -Zbin(Xs,t,T,W,index,1,bpN,bN);
            cons = [cons, Wc'* Lamda_s1 == Zs1'];
            cons = [cons, Uc'* Lamda_s1 + Lamda2_s1== 0];
            cons = [cons, hc'* Lamda_s1 >= 0];
            cons = [cons, Lamda2_s1 >= 0];
            %d>0 constraints
            Lamda_d0 = [sdpvar(t*W*kl + 3*t,1); sdpvar(t*W + 3*t,1)];
            Lamda2_d0 = sdpvar((t*W)*vN + 3*t,1);
            Zd0 = Zbin(Xd,t,T,W,index,0,bpN,bN);
            cons = [cons, Wc'* Lamda_d0 == Zd0'];
            cons = [cons, Uc'* Lamda_d0 + Lamda2_d0== 0];
            cons = [cons, hc'* Lamda_d0 >= 0];
            cons = [cons, Lamda2_d0 >= 0];
            %d<1
            Lamda_d1 = [sdpvar(t*W*kl + 3*t,1); sdpvar(t*W + 3*t,1)];
            Lamda2_d1 = sdpvar((t*W)*vN + 3*t,1);
            Zd1 = -Zbin(Xd,t,T,W,index,1,bpN,bN);
            cons = [cons, Wc'* Lamda_d1 == Zd1'];
            cons = [cons, Uc'* Lamda_d1 + Lamda2_d1== 0];
            cons = [cons, hc'* Lamda_d1 >= 0];
            cons = [cons, Lamda2_d1 >= 0];
            omiga = max(0,t-Ton(index))+1;
            if (omiga<=t)
                Lamda_on = [sdpvar(t*W*kl + 3*t,1); sdpvar(t*W + 3*t,1)];
                Lamda2_on = sdpvar((t*W)*vN + 3*t,1);
                Bon_s = zeros(1,1+t*W*bpN);
                for o = omiga:t
                    Bon_s = Bon_s + [Xs(index,o),Xs(index,T+1+Sumt(o-1)*W*bpN:T+Sumt(o)*W*bpN), zeros(1,(t-o)*W*bpN)];
                end
                Bon_u = [Xu(index,t),Xu(index,T+1+Sumt(t-1)*W*bpN:T+Sumt(t)*W*bpN)];
                Zmid = Bon_u - Bon_s; 
                Zon = [];
                for tt = 1:t
                    if tt == t
                        Zon = [Zon,0,0,Zmid(1,1)];
                    else
                        Zon = [Zon,0,0,0];
                    end
                end
                for tt = 1:t 
                    for w = 1:W
                        Zon = [Zon,0,zeros(1,bN),Zmid(1,2+(tt-1)*W*bpN+(w-1)*bpN:1+(tt-1)*W*bpN+w*bpN)];
                    end
                end
                cons = [cons, Wc'* Lamda_on == Zon'];
                cons = [cons, Uc'* Lamda_on + Lamda2_on== 0];
                cons = [cons, hc'* Lamda_on >= 0];
                cons = [cons, Lamda2_on >= 0];
            end
            %shutdown constraints
            omiga = max(0,t-Toff(index))+1;
            if (t >= 1+L0(index)&&(omiga<=t)) 
                Lamda_off = [sdpvar(t*W*kl + 3*t,1); sdpvar(t*W + 3*t,1)];
                Lamda2_off = sdpvar((t*W)*vN + 3*t,1);
                Boff_d = zeros(1,1+t*W*bpN);
                for o = omiga:t
                    Boff_d = Boff_d - [Xd(index,o), Xd(index,T+1+Sumt(o-1)*W*bpN:T+Sumt(o)*W*bpN), zeros(1,(t-o)*W*bpN)];
                end
                Boff_u = -[Xu(index,t),Xu(index,T+1+Sumt(t-1)*W*bpN:T+Sumt(t)*W*bpN)];
                Boff = Boff_u + Boff_d;
                Zoff = [];
                for tt = 1:t
                    if tt == t
                        Zoff = [Zoff,0,0,1+Boff(1,1)];
                    else
                        Zoff = [Zoff,0,0,0];
                    end
                end
                for tt = 1:t 
                    for w = 1:W
                        Zoff = [Zoff,0,zeros(1,bN),Boff(1,2+(tt-1)*W*bpN+(w-1)*bpN:1+(tt-1)*W*bpN+w*bpN)];
                    end
                end
                cons = [cons, Wc'* Lamda_off == Zoff'];
                cons = [cons, Uc'* Lamda_off + Lamda2_off== 0];
                cons = [cons, hc'* Lamda_off >= 0];
                cons = [cons, Lamda2_off >= 0]; 
            end
            %end of only binary
        end
%         % transmission limitation
        for i = 1:size(all_branch.I,1)
            left = all_branch.I(i);
            right = all_branch.J(i);
            abs_x4branch = abs(1/B(left,right));  % the absolute value of the impedance of the current branch |x_ij|
            % right side of the inequility
            Lamda_Fup = [sdpvar(t*W*kl + 3*t,1); sdpvar(t*W + 3*t,1)];
            Lamda2_Fup = sdpvar((t*W)*vN + 3*t,1);
            AFup = [Ytheta(right,t) - Ytheta(left,t) + abs_x4branch*all_branch.P(i),Ytheta(right,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN) - Ytheta(left,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN)]; % theta_j - theta_i + F
            ZFup = [];
            for tt = 1:t
                if tt == t
                    ZFup = [ZFup,0,AFup(1,1),0];
                else
                    ZFup = [ZFup,0,0,0];
                end
            end
            for tt = 1:t
                for w = 1:W
                    ZFup = [ZFup,0,AFup(1,2+(tt-1)*W*bN+(w-1)*bN:1+(tt-1)*W*bN+w*bN),zeros(1,bpN)];
                end
            end
            cons = [cons, Wc'* Lamda_Fup == ZFup'];
            cons = [cons, Uc'* Lamda_Fup + Lamda2_Fup == 0];
            cons = [cons, hc'* Lamda_Fup >= 0];
            cons = [cons, Lamda2_Fup >= 0];
            % left side
            Lamda_Flow = [sdpvar(t*W*kl + 3*t,1); sdpvar(t*W + 3*t,1)];
            Lamda2_Flow = sdpvar((t*W)*vN + 3*t,1);
            AFlow = [Ytheta(left,t) - Ytheta(right,t) + abs_x4branch*all_branch.P(i),Ytheta(left,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN) - Ytheta(right,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN)]; % theta_i - theta_j + F
            ZFlow = [];
            for tt = 1:t
                if tt == t
                    ZFlow = [ZFlow,0,AFlow(1,1),0];
                else
                    ZFlow = [ZFlow,0,0,0];
                end
            end
            for tt = 1:t
                for w = 1:W
                    ZFlow = [ZFlow,0,AFlow(1,2+(tt-1)*W*bN+(w-1)*bN:1+(tt-1)*W*bN+w*bN),zeros(1,bpN)];
                end
            end
            cons = [cons, Wc'* Lamda_Flow == ZFlow'];
            cons = [cons, Uc'* Lamda_Flow + Lamda2_Flow == 0];
            cons = [cons, hc'* Lamda_Flow >= 0];
            cons = [cons, Lamda2_Flow >= 0];
        end
    end
end
% end of construct inequality constraints

% construct equality constraints
% on/off unit state
for g = 1:G
    for t = 1:T
        if t == 1
            Zmot_0 = [0,0,Xs(g,t)] - [0,0,Xd(g,t)] - [0,0,Xu(g,t)];
            Zmot = [0,zeros(1,bN),Xs(g,T+1:T+W*bpN)] - [0,zeros(1,bN),Xd(g,T+1:T+W*bpN)] - [0,zeros(1,bN),Xu(g,T+1:T+W*bpN)];
        else
            Zmot_0 = [0,0,Xs(g,t)-Xd(g,t)-Xu(g,t)+Xu(g,t-1)];
            Bmot = Xs(g,T+1+Sumt(t-1)*W*bpN:T+Sumt(t)*W*bpN) - Xd(g,T+1+Sumt(t-1)*W*bpN:T+Sumt(t)*W*bpN) - Xu(g,T+1+Sumt(t-1)*W*bpN:T+Sumt(t)*W*bpN) + [Xu(g,T+1+Sumt(t-2)*W*bpN:T+Sumt(t-1)*W*bpN),zeros(1,W*bpN)]; 
            Zmot = [];
            for tt = 1:t
                for w = 1:W
                    Zmot = [Zmot,0,zeros(1,bN),Bmot(1,1+(tt-1)*W*bpN+(w-1)*bpN:(tt-1)*W*bpN+w*bpN)];
                end
            end
        end
        cons = [cons, Zmot_0*ones(3,1) == 0]; 
        cons = [cons, Zmot == 0]; 
    end
end

% DC network constraints
for n = 1:N 
    for t = 1:T
        Zac = zeros(1,1+t*W*bN); %continuous decision variable
        dac = 0; %constant term
        Zac = Zac + [Yl(n,t),Yl(n,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN)];
        if ismember(n,Gbus) %if thermal units
            index = find(Gbus==n);
            Zac = Zac + [YG(index,t),YG(index,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN)];
        end
        if ismember(n,Wbus) %if wind units
            index = find(Wbus==n);
            Zac = Zac + [YPw(index,t),YPw(index,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN)];
        end
        for j = 1:N 
            if B(n,j) ~=0
                Zac = Zac - B(n,j) * [Ytheta(j,t),Ytheta(j,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN)]; 
            end
        end
        if ismember(n,Dbus) %if load bus
            index = find(Dbus==n);
            dac = dac + SCUC_data.busLoad.node_P(t,index);
        end
        cons = [cons, Zac(1,1) == dac]; 
        cons = [cons, Zac(2:t*W*bN)' == 0]; 
    end
end
% end of DC network constraints
% end of construct equality constraints

% construct object
Cobj = zeros(1,T*W); %ksi coefficient matrix
Aobj = zeros(1,T*W*bN); %ksi_hat coefficient matrix
Bobj = zeros(1,T*W*bpN); %ksi_tine coefficient matrix
Cobj0 = zeros(1,T); %ksi0 coefficient matrix
Aobj0 = zeros(1,T); %ksi_hat0 coefficient matrix
Bobj0 = zeros(1,T); %ksi_tine0 coefficient matrix
Uobj_l = []; %convex hull's verts matrix of objective
Uobj_c = [];
for t = 1:T
    for n = 1:N
        Aobj0 = Aobj0+ 100*[zeros(1,t-1),Yl(n,t),zeros(1,T-t)];
        Aobj = Aobj + 100*[Yl(n,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN),zeros(1,(T-t)*W*bN)]; 
        if ismember(n,Gbus) 
            index = find(Gbus==n);
            Aobj0 = Aobj0+ [zeros(1,t-1),Yz(index,t),zeros(1,T-t)] + [zeros(1,t-1),YS(index,t),zeros(1,T-t)];
            Bobj0 = Bobj0 + [zeros(1,t-1),Xs(index,t),zeros(1,T-t)];
            Aobj = Aobj + [Yz(index,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN) + YS(index,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN), zeros(1,(T-t)*W*bN)]; 
            Bobj = Bobj + [Chot(index) * Xs(index,1+Sumt(t-1)*W*bpN:Sumt(t)*W*bpN), zeros(1,(T-t)*W*bpN)]; 
        end
    end 
    Ulamda = blkdiag(verts(:,1+(t-1)*2*vN:vN+(t-1)*2*vN),verts(:,1+vN+(t-1)*2*vN:2*vN*t));
    Uconstant = []; 
    for w = 1:W
        Uconstant = blkdiag(Uconstant,ones(1,vN));
    end
    Uobj_l = blkdiag(Ulamda_0,Uobj_l);
    Uobj_l = blkdiag(Uobj_l,-Ulamda);
    Uobj_c = blkdiag(Uconstant_0,Uobj_c);
    Uobj_c = blkdiag(Uobj_c,Uconstant);
end

U = [Uobj_l;Uobj_c];
W_l = diag(ones(T*W*kl+3*T,1));
W_c = zeros(T*W+3*T,T*W*kl+3*T);
W_obj = [W_l;W_c];
Zobj = [];
Zobj_0 = [];
for t = 1:T
    Zobj_0 = [Zobj_0,0,Aobj0(1,t),Bobj0(1,t)];
    for w = 1:W
        Zobj = [Zobj,0,Aobj(1,1+(t-1)*2*bN+(w-1)*bN:(t-1)*2*bN+w*bN),Bobj(1,1+(t-1)*2*bpN+(w-1)*bpN:(t-1)*2*bpN+w*bpN)];
    end
end
Zobj = [Zobj_0,Zobj]; %ksi0,ksi_hat0,ksi_tine0,ksi,ksi_hat,ksi_tine
hobj = [zeros(1,T*W*kl+3*T),ones(1,T*W+3*T)]'; 
v = sdpvar(1,M); %first lagrangian multipliers v
beta = sdpvar(1); %second lagrangian multipliers beta
y1 = sdpvar((T*W*(kl+1))+2*3*T,M);
constraints = [];
constraints = [constraints, beta >= 0]; 
for m = 1:M 
    y_j = y1(:,m); 
    mid_obj = Zobj' - W_obj'*y_j;
    constraints = [constraints, (mid_obj' * ksi_wan{m} + hobj' * y_j) <= M*v(1,m)];
    constraints = [constraints, abs(mid_obj) <= beta];
    constraints = [constraints, U' * y_j >= 0];
end
% end of construct object

constraints = [constraints, cons];

objective = sum(v) + eps*beta; %setting objective

% options = sdpsettings('verbose',2,'debug',1,'savesolveroutput',1,'savesolverinput',1);
% options = sdpsettings('verbose',2,'solver','mosek','debug',1,'savesolveroutput',1,'savesolverinput',1);
% options.mosek.MSK_DPAR_OPTIMIZER_MAX_TIME = 1000;
options = sdpsettings('verbose',2,'solver','cplex','debug',1,'savesolveroutput',1,'savesolverinput',1);
% options.cplex.timelimit = 1000;
sol = optimize(constraints,objective,options);

% Analyze error flags
if sol.problem == 0
    % Extract and display value
    Obj = value(objective);
    disp(Obj);

    real_PG = zeros(G,T); %thermal generation
    real_z = zeros(G,T); %generation cost of thermal units
    real_S = zeros(G,T); %startpu cost of thermal units
    real_l = zeros(N,T); %load loss
    real_PW = zeros(W,T); %wind generation
    real_theta = zeros(N,T); %phase angle
    real_u = zeros(G,T); %units state
    real_s = zeros(G,T); %turn on
    real_d = zeros(G,T); %turn off
   
    ksi_tine = ksi_wan2{18}(1+T*W*bN:T*W*bN+T*W*bpN);
    ksi_hat = ksi_wan2{18}(1:T*W*bN);
    for g = 1:G
        for t = 1:T
            real_PG(g,t) = YG(g,t) + sum(value(YG(g,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN))' .* ksi_hat(1:t*W*bN));
            real_z(g,t) = Yz(g,t) + sum(value(Yz(g,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN))' .* ksi_hat(1:t*W*bN));
            real_S(g,t) = YS(g,t) + sum(value(YS(g,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN))' .* ksi_hat(1:t*W*bN));
            real_u(g,t) = Xu(g,t) + sum(value(Xu(g,T+1+Sumt(t-1)*W*bpN:T+Sumt(t)*W*bpN))' .* ksi_tine(1:t*W*bpN));
            real_s(g,t) = Xs(g,t) + sum(value(Xs(g,T+1+Sumt(t-1)*W*bpN:T+Sumt(t)*W*bpN))' .* ksi_tine(1:t*W*bpN));
            real_d(g,t) = Xd(g,t) + sum(value(Xd(g,T+1+Sumt(t-1)*W*bpN:T+Sumt(t)*W*bpN))' .* ksi_tine(1:t*W*bpN));
        end
    end
    
    for w = 1:W
        for t = 1:T
            real_PW(w,t) = YPw(w,t) + sum(value(YPw(w,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN))' .* ksi_hat(1:t*W*bN));
        end
    end

    for n = 1:N
        for t = 1:T
            real_theta(n,t) = Ytheta(n,t) + sum(value(Ytheta(n,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN))' .* ksi_hat(1:t*W*bN));
            real_l(n,t) = Yl(n,t) + sum(value(Yl(n,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN))' .* ksi_hat(1:t*W*bN));
        end
    end
    
    real_yyy = zeros(N,T);
    for n = 1:N 
        for t = 1:T
            real_yyy(n,t) = real_yyy(n,t) + real_l(n,t);
            if ismember(n,Gbus) 
                index = find(Gbus==n);
                real_yyy(n,t) = real_yyy(n,t) + real_PG(index,t);
            end
            if ismember(n,Wbus)
                index = find(Wbus==n);
                real_yyy(n,t) = real_yyy(n,t) + real_PW(index,t);
            end
            for j = 1:N 
                if B(n,j) ~=0
                    real_yyy(n,t) = real_yyy(n,t) - B(n,j)*real_theta(j,t);
                end
            end
            if ismember(n,Dbus)
                index = find(Dbus==n);
                real_yyy(n,t) = real_yyy(n,t) - SCUC_data.busLoad.node_P(t,index);
            end
        end
    end
    disp(sol.solvertime);
else
    disp('Oh shit!, something was wrong!');
    sol.info
    yalmiperror(sol.problem)
end
