%%  BDR_UCC_DRO_Last  by yy 2022.04.10

FileName_SCUC = 'data/SCUC6.txt';
FileName_Wind = 'data\Wind_power.xlsx';

SCUC_data = ReadDataSCUC(FileName_SCUC);
Wind_data = ReadWindData(SCUC_data,FileName_Wind);

% T = SCUC_data.totalLoad.T;  % 时段数T
T = 12;  % 时段数T
G = SCUC_data.units.N;      % 火力发电机数
N = SCUC_data.baseparameters.busN;  % 节点总数
W = SCUC_data.Wind.wind_Number; % 风电机组个数
M = size(Wind_data{1,1},2); % 样本数量
Gbus = SCUC_data.units.bus_G; % 火力发电机节点号
Wbus = SCUC_data.Wind.wind_Node; % 风电机组节点号
Dbus = SCUC_data.busLoad.bus_PDQR; % 负载节点号
Pramp = SCUC_data.units.ramp; %爬坡约束
Pstart = 0.5 * SCUC_data.units.PG_up; %开机爬约束
Pshut = SCUC_data.units.PG_up; %关机爬约束 
Pup = SCUC_data.units.PG_up; %出力上界
Plow = SCUC_data.units.PG_low; %出力下界
Ton = SCUC_data.units.T_on; %最小开机时间
Toff = SCUC_data.units.T_off; %最小关机时间
L = 4; %给定的参数，用于线性化发电费用
fa = SCUC_data.units.alpha; %发电费用常数项
fb = SCUC_data.units.beta; %发电费用一次项
fc = SCUC_data.units.gamma; %发电费用二次项
Ccold = SCUC_data.units.start_cost*2; %冷启动费用
Chot = SCUC_data.units.start_cost; % 热启动费用
Tcold = min(max(Toff,1),T); %冷启动时间
Lup = 0.3*ones(N,1); % 功率切割上界
Llow = 0*ones(N,1); %功率切割下界

% 形成直流潮流系数矩阵B
type_of_pf = 'DC';
Y = SCUC_nodeY(SCUC_data,type_of_pf);
B = -Y.B; %因为是直流方程 所以B忽略了电阻 只考虑电抗

u0 = zeros(G,1); %机组初始状态，目前假设都是关机
T0 = zeros(G,1); %机组初始时还需开机或关机时段数
P0 = zeros(G,1); %对应初始发电也为0
U0 = zeros(G,1); 
L0 = zeros(G,1); 
tao0 = Toff+Tcold+1; %冷启动需要时间？

for index = 1:G
    mid_u = min(T,u0(index) * (Ton(index) - T0(index)));
    mid_l = min(T,(1-u0(index)) * (Toff(index) + T0(index)));
    U0(index) = max(0,mid_u);
    L0(index) = max(0,mid_l);
end

all_branch.I = [ SCUC_data.branch.I; SCUC_data.branchTransformer.I ]; % 所有支路起点 前是支路起点 后是变压器支路起点
all_branch.J = [ SCUC_data.branch.J; SCUC_data.branchTransformer.J ]; % 所有支路终点
all_branch.P = [ SCUC_data.branch.P; SCUC_data.branchTransformer.P ]; % 支路功率上限

% 计算 uncertainty set
PW = []; % 将风力发电整合到一个大矩阵中，每列按先变t后i排列
for i = 1:W
    PW = [PW; Wind_data{1,i}];
end
PW = [PW(1:T,:);PW(24+1:24+T,:)];
% for i = 1:W
%     PW = [PW; 0.05*ones(24,30)];
% end
% 风力发电均值
PW_miu = zeros(W*T,1);
for i = 1:M
    PW_miu = PW_miu + PW(:,i);
end
PW_miu = PW_miu / M;
% 风电偏差
VW = zeros(W*T,M);
for j = 1:M
    VW(:,j) =  PW(:,j) - PW_miu;
end
VW_positive =  max(VW,[],2);
VW_negative =  min(VW,[],2);
% ksi 支撑集
ksi_positive = [1;VW_positive];
ksi_negative = [1;VW_negative];
% end of 计算 uncertainty set

% non-linear decision rule
bN = 3; % 分段数
bpN = bN-1; % 断点数
bp = zeros(bN - 1,1); % 等分点向量
kl = 1+bpN+bN; %一个ksi的维度
verts = []; % set of vertices
for t = 1:T
    for i = 1:W
        step = (VW_positive((t-1)*W+i) - VW_negative((t-1)*W+i)) / bN; %计算每段区间长度
        bp = linspace(VW_negative((t-1)*W+i)+ step,VW_positive((t-1)*W+i)-step,bN-1)'; %从[a+s,b-s]中间区间取等分点，因为linspace函数只处理闭区间
        verts = [verts, vers(bN,bp,VW_positive((t-1)*W+i),VW_negative((t-1)*W+i))];
    end
end
verts_nonksi = verts(2:6,:);
vN = length(verts)/(W*T); %结点数
% verts = [ones(2*bN,(bN+1)*W), verts];
% end of non-linear decision rule

% 将所有ksi做决策处理
for m = 1:M
    ksi = VW(:,m);
    ksi_hat = [];
    ksi_tine = [];
    ksi_hat2 = [];
    ksi_tine2 = [];
    ksi_wan_now = [];
    for t = 1:T
        for i = 1:W
            step = (VW_positive((t-1)*W+i) - VW_negative((t-1)*W+i)) / bN; %计算每段区间长度
            bp = linspace(VW_negative((t-1)*W+i)+ step,VW_positive((t-1)*W+i)-step,bN-1)'; %从[a+s,b-s]中间区间取等分点，因为linspace函数只处理闭区间
            bp = [VW_negative((t-1)*W+i);bp;VW_positive((t-1)*W+i)];
            ksi_hat = L_hat(ksi((t-1)*W+i),bN,bp);
            ksi_tine = L_tine(ksi((t-1)*W+i),bN,bp);
            ksi_hat2 = [ksi_hat2; L_hat(ksi((i-1)*T+t),bN,bp)];
            ksi_tine2 = [ksi_tine2;  L_tine(ksi((i-1)*T+t),bN,bp)];
            ksi_wan_now = [ksi_wan_now;ksi((t-1)*W+i);ksi_hat;ksi_tine];
        end
    end
%     ksi_wan{m} = [ones(T,1); ones(T,1); ksi_hat; ksi_tine]; %按ksi_hat0;ksi_tine0;ksi_hat;ksi_tine排列
    ksi_wan{m} = [ ones(T,1);ones(T,1);ones(T,1);ksi_wan_now]; %按ksi0;ksi_hat0;ksi_tine0;ksi;ksi_hat;ksi_tine排列
    ksi_wan2{m} = [ ksi_hat2; ksi_tine2]; %按ksi_hat;ksi_tine排列
end
% end of 将所有ksi做决策处理

% % 构造wasserstein模糊集
% eta = 0.95; %设置置信度
% rho = sdpvar(1); %模糊集半径
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
% % end of 构造wasserstein模糊集
eps = 4.5;
% 构造闭凸包
% 极点在上节求好了 verts

% end of 构造闭凸包

%计算ksi长度总长度
leng_k = 0;
for k = 1:T
    leng_k = leng_k+k;
end
%end of 计算ksi长度


% 构造约束
% 连续变量定义 每个时段变量数不同，在用到时定义
YG = sdpvar(G,T+W*bN*leng_k); %火力机组出力决策变量，系数矩阵按Y0,分段，风电机组，时段排列
% YG = 1:64;
% YG = [YG;YG;YG;YG;YG;YG];
Ytheta = sdpvar(N,T+W*bN*leng_k); %节点相角决策变量，每个节点都有theta，不止发电机节点有，所以有N行
Yz = sdpvar(G,T+W*bN*leng_k); % 火力机组出力费用决策变量
YS = sdpvar(G,T+W*bN*leng_k); % 火力机组开机费用决策变量
YPw = sdpvar(W,T+W*bN*leng_k); %W行，W个风电机组
Yl = sdpvar(N,T+W*bN*leng_k); %切负载决策变量
% 整型变量
% 是取值为-1，0，1的整型变量，不是0，1binary变量
Xu = intvar(G,T+W*bpN*leng_k); %机组状态决策变量
Xs = intvar(G,T+W*bpN*leng_k); %开机动作决策变量
Xd = intvar(G,T+W*bpN*leng_k); %关机动作决策变量

% 开始约束构造
cons = [];
cons = [cons, -1 <= Xu <= 1]; % 将整型决策变量约束在-1，0，1
cons = [cons, -1 <= Xs <= 1];
cons = [cons, -1 <= Xd <= 1];

% 参考节点
%每个模型都只设定第一个节点为参考节点？
Ytheta(1,:) = zeros(1,T+W*bN*leng_k); %每个节点都有theta，不止发电机节点有，所有有N行
% end of 参考节点

% 初始状态保持
for g = 1:G % 设置了开始都是关机，所以这里保持关机状态。没考虑初始是开机情况
    for t = 1:(U0(g)+L0(g))
        Xu(g,T+1+Sumt(t-1)*W*bpN:T+Sumt(t)*W*bpN) = zeros(1,t*W*bpN); %将需要保持关机态时段Xu的决策变量系数置零
        Xu(g,t) = 0; %Xu0 也要置零
    end
end
% end of 初始状态保持

% 不等式约束统一构造
% 构造通用的U矩阵，h向量 公用常量
Uc_l = []; %迭代生成顶点矩阵
Uc_c = []; %lamda系数矩阵
Uc = [];
Uconstant_0 = []; %顶点系数lamda组合为1约束的系数矩阵 for ksi0
Uconstant_0 = blkdiag(1, Uconstant_0);%在最开始为Uconstant矩阵补上系数1 for ksi0 常数项
Uconstant_0 = blkdiag(1, Uconstant_0);%在最开始为Uconstant矩阵补上系数1 for ksi_hat0 常数项
Uconstant_0 = blkdiag(1, Uconstant_0);%在最开始为Uconstant矩阵补上系数1 for ksi_tine0 常数项
Ulamda_0 = [];
Ulamda_0 = blkdiag(-1, Ulamda_0); %在最开始为Ulamda矩阵补上系数-1 for ksi0 常数项
Ulamda_0 = blkdiag(-1, Ulamda_0); %在最开始为Ulamda矩阵补上系数-1 for ksi_hat0 常数项
Ulamda_0 = blkdiag(-1, Ulamda_0); %在最开始为Ulamda矩阵补上系数-1 for ksi_tine0 常数项
for t = 1:T
    Uconstant = []; %顶点系数lamda组合为1约束的系数矩阵
    for w = 1:W
        Uconstant = blkdiag(Uconstant,ones(1,vN));
    end
    Ulamda = blkdiag(-verts(:,1+(t-1)*2*vN:vN+(t-1)*2*vN),-verts(:,1+vN+(t-1)*2*vN:2*vN*t)); %顶点凸组合表示ksi约束的系数矩阵
    Uc_l = blkdiag(Ulamda_0,Uc_l); % U矩阵上部分，由凸包顶点构成
    Uc_l = blkdiag(Uc_l,Ulamda); % U矩阵上部分，由凸包顶点构成
    Uc_c = blkdiag(Uconstant_0,Uc_c); %U矩阵下部分，由1组成实现顶点系数求和为1
    Uc_c = blkdiag(Uc_c,Uconstant); %U矩阵下部分，由1组成实现顶点系数求和为1
    Uc = [Uc_l;Uc_c]; %凸包等式约束U矩阵
    Wc_l = diag(ones(t*W*kl+3*t,1)); %W矩阵上部分，单位阵
    Wc_c = zeros(t*W+3*t,t*W*kl+3*t); %W矩阵下部分，0矩阵
    Wc = [Wc_l;Wc_c];
    hc = [zeros(t*W*kl + 3*t,1);ones(t*W+3*t,1)]; %构造常数项向量 = 两倍的Z向量长度
    for n = 1:N %遍历所有节点
        % 功率切割上界约束
        Lamda_lup = [sdpvar(t*W*kl + 3*t,1); sdpvar(t*W + 3*t,1)]; %构造等式对偶乘子lamda向量 = h向量维度
        Lamda2_lup = sdpvar((t*W)*vN + 3*t,1); %构造不等式对偶乘子lamda向量 = verts维度
        Zlup = -Zhat(Yl,t,T,W,n,Lup(n,1),bpN,bN);
        cons = [cons, Wc'* Lamda_lup == Zlup'];
        cons = [cons, Uc'* Lamda_lup + Lamda2_lup == 0];
        cons = [cons, hc'* Lamda_lup >= 0];
        cons = [cons, Lamda2_lup >= 0];
        % 功率切割下界约束
        Lamda_llow = [sdpvar(t*W*kl + 3*t,1); sdpvar(t*W + 3*t,1)]; %构造等式对偶乘子lamda向量 = h向量维度
        Lamda2_llow = sdpvar((t*W)*vN + 3*t,1); %构造不等式对偶乘子lamda向量 = verts维度
        Zllow = Zhat(Yl,t,T,W,n,Llow(n,1),bpN,bN);
        cons = [cons, Wc'* Lamda_llow == Zllow'];
        cons = [cons, Uc'* Lamda_llow + Lamda2_llow == 0];
        cons = [cons, hc'* Lamda_llow >= 0];
        cons = [cons, Lamda2_llow >= 0];
        %判断是否是风电发力节点，加上风电发力约束
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
            %风电出力下界
            Lamda_wlow = [sdpvar(t*W*(1+bN+bpN) + 3*t,1); sdpvar(t*W + 3*t,1)]; %构造等式对偶乘子lamda向量 = h向量维度
            Lamda2_wlow = sdpvar((t*W)*vN + 3*t,1); %构造不等式对偶乘子lamda向量 = verts维度
            Zwlow = Zhat(YPw,t,T,W,index,0,bpN,bN);
            cons = [cons, Wc'* Lamda_wlow == Zwlow'];
            cons = [cons, Uc'* Lamda_wlow + Lamda2_wlow == 0];
            cons = [cons, hc'* Lamda_wlow >= 0];
            cons = [cons, Lamda2_wlow >= 0];
        end
        %判断是否是火力发电节点，是的话，加上相关约束
        if ismember(n,Gbus) 
            index = find(Gbus==n);
            %最低发电费用约束
            Lamda_wan = [sdpvar(t*W*kl + 3*t,1); sdpvar(t*W + 3*t,1)];
            Lamda2_wan = sdpvar((t*W)*vN + 3*t,1);
            Zwan = Zhat(YS,t,T,W,index,0,bpN,bN);
            cons = [cons, Wc'* Lamda_wan == Zwan'];
            cons = [cons, Uc'* Lamda_wan + Lamda2_wan == 0];
            cons = [cons, hc'* Lamda_wan >= 0];
            cons = [cons, Lamda2_wan >= 0];
            
            % 发电费用线性化
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
            
            % 开机费用约束
            if (t-tao0(index)+1<=0)&&(max(0,-T0(index))<abs(t-tao0(index)+1))
                fit = 1;
            else
                fit = 0;
            end
            tao_it = max(1,t-Toff(index)-Tcold(index));
            Lamda_open = [sdpvar(t*W*kl + 3*t,1); sdpvar(t*W + 3*t,1)];
            Lamda2_open = sdpvar((t*W)*vN + 3*t,1);
            Ccost = Chot(index) - Ccold(index); %计算费用公共系数
            Aopen = [YS(index,t)-fit*Ccost,YS(index,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN)];
            Bopen = zeros(1,1+t*W*bpN);
            for y = tao_it : t-1 % 构造关机状态dit的决策变量
                Bopen = Bopen - [Xd(index,y),Xd(index,T+1+Sumt(y-1)*W*bpN:T+Sumt(y)*W*bpN), zeros(1,(t-y)*W*bpN)];
            end
            Bopen = Bopen + [Xs(index,t),Xs(index,T+1+Sumt(t-1)*W*bpN:T+Sumt(t)*W*bpN)]; %加上开机状态sit的决策变量
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
            cons = [cons, Lamda2_open >= 0]; %不等式constraints乘子有约束，等式的无约束
            
            %出力上界
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
            for tt = 1:t %按当前时段用到的前t个时段W个机组取出对应系数
                for w = 1:W
                    Zup = [Zup,0,-Aup(1,2+(tt-1)*W*bN+(w-1)*bN:1+(tt-1)*W*bN+w*bN),Bup(1,2+(tt-1)*W*bpN+(w-1)*bpN:1+(tt-1)*W*bpN+w*bpN)];
                end
            end
            cons = [cons, Wc'* Lamda_up == Zup'];
            cons = [cons, Uc'* Lamda_up + Lamda2_up == 0];
            cons = [cons, hc'* Lamda_up >= 0];
            cons = [cons, Lamda2_up >= 0]; 
            
            %出力下界
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
            
            %上爬坡
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
            
            %下爬坡
            Lamda_ramp_down = [sdpvar(t*W*kl + 3*t,1); sdpvar(t*W + 3*t,1)];
            Lamda2_ramp_down = sdpvar((t*W)*vN + 3*t,1);
            if t == 1
                Aramp_down = [YG(index,t),YG(index,T+1:T+W*bN)]; %当t=1时，剩下d,u和p1三项
                Bramp_down = Pshut(index) * [Xd(index,t),Xd(index,T+1:T+W*bpN)] + Pramp(index) * [Xu(index,t),Xu(index,T+1:T+W*bpN)];
            else
                Aramp_down = [YG(index,t),YG(index,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN)] - [YG(index,t-1),YG(index,T+1+Sumt(t-2)*W*bN:T+Sumt(t-1)*W*bN), zeros(1,W*bN)]; %构造连续变量ksi_hat前系数 Pt - Pt-1
                Bramp_down = Pshut(index) * [Xd(index,t),Xd(index,T+1+Sumt(t-1)*W*bpN:T+Sumt(t)*W*bpN)] + Pramp(index) * [Xu(index,t),Xu(index,T+1+Sumt(t-1)*W*bpN:T+Sumt(t)*W*bpN)]; %构造整型ksi_tine前系数
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
            %开机状态u>0约束
            Lamda_u0 = [sdpvar(t*W*kl + 3*t,1); sdpvar(t*W + 3*t,1)];
            Lamda2_u0 = sdpvar((t*W)*vN + 3*t,1);
            Zu0 = Zbin(Xu,t,T,W,index,0,bpN,bN);
            cons = [cons, Wc'* Lamda_u0 == Zu0'];
            cons = [cons, Uc'* Lamda_u0 + Lamda2_u0== 0];
            cons = [cons, hc'* Lamda_u0 >= 0];
            cons = [cons, Lamda2_u0 >= 0];
            %u<1限制
            Lamda_u1 = [sdpvar(t*W*kl + 3*t,1); sdpvar(t*W + 3*t,1)];
            Lamda2_u1 = sdpvar((t*W)*vN + 3*t,1);
            Zu1 = -Zbin(Xu,t,T,W,index,1,bpN,bN);
            cons = [cons, Wc'* Lamda_u1 == Zu1'];
            cons = [cons, Uc'* Lamda_u1 + Lamda2_u1== 0];
            cons = [cons, hc'* Lamda_u1 >= 0];
            cons = [cons, Lamda2_u1 >= 0];
            %开机动作s>0约束
            Lamda_s0 = [sdpvar(t*W*kl + 3*t,1); sdpvar(t*W + 3*t,1)];
            Lamda2_s0 = sdpvar((t*W)*vN + 3*t,1);
            Zs0 = Zbin(Xs,t,T,W,index,0,bpN,bN);
            cons = [cons, Wc'* Lamda_s0 == Zs0'];
            cons = [cons, Uc'* Lamda_s0 + Lamda2_s0== 0];
            cons = [cons, hc'* Lamda_s0 >= 0];
            cons = [cons, Lamda2_s0 >= 0];
            %s<1限制
            Lamda_s1 = [sdpvar(t*W*kl + 3*t,1); sdpvar(t*W + 3*t,1)];
            Lamda2_s1 = sdpvar((t*W)*vN + 3*t,1);
            Zs1 = -Zbin(Xs,t,T,W,index,1,bpN,bN);
            cons = [cons, Wc'* Lamda_s1 == Zs1'];
            cons = [cons, Uc'* Lamda_s1 + Lamda2_s1== 0];
            cons = [cons, hc'* Lamda_s1 >= 0];
            cons = [cons, Lamda2_s1 >= 0];
            %关机动作d>0约束
            Lamda_d0 = [sdpvar(t*W*kl + 3*t,1); sdpvar(t*W + 3*t,1)];
            Lamda2_d0 = sdpvar((t*W)*vN + 3*t,1);
            Zd0 = Zbin(Xd,t,T,W,index,0,bpN,bN);
            cons = [cons, Wc'* Lamda_d0 == Zd0'];
            cons = [cons, Uc'* Lamda_d0 + Lamda2_d0== 0];
            cons = [cons, hc'* Lamda_d0 >= 0];
            cons = [cons, Lamda2_d0 >= 0];
            %d<1限制
            Lamda_d1 = [sdpvar(t*W*kl + 3*t,1); sdpvar(t*W + 3*t,1)];
            Lamda2_d1 = sdpvar((t*W)*vN + 3*t,1);
            Zd1 = -Zbin(Xd,t,T,W,index,1,bpN,bN);
            cons = [cons, Wc'* Lamda_d1 == Zd1'];
            cons = [cons, Uc'* Lamda_d1 + Lamda2_d1== 0];
            cons = [cons, hc'* Lamda_d1 >= 0];
            cons = [cons, Lamda2_d1 >= 0];
            %开机限制
            %由于设置了U0都为0，所以这版代码没有对t进行像关机限制时那样的判断
            omiga = max(0,t-Ton(index))+1;
            if (omiga<=t)
                Lamda_on = [sdpvar(t*W*kl + 3*t,1); sdpvar(t*W + 3*t,1)];
                Lamda2_on = sdpvar((t*W)*vN + 3*t,1);
                Bon_s = zeros(1,1+t*W*bpN);
                for o = omiga:t
                    Bon_s = Bon_s + [Xs(index,o),Xs(index,T+1+Sumt(o-1)*W*bpN:T+Sumt(o)*W*bpN), zeros(1,(t-o)*W*bpN)];
                end
                Bon_u = [Xu(index,t),Xu(index,T+1+Sumt(t-1)*W*bpN:T+Sumt(t)*W*bpN)];
                Zmid = Bon_u - Bon_s; %构造当前时段当前机组用到的u-s向量
                Zon = [];
                for tt = 1:t
                    if tt == t
                        Zon = [Zon,0,0,Zmid(1,1)];
                    else
                        Zon = [Zon,0,0,0];
                    end
                end
                for tt = 1:t %按当前时段用到的前t个时段W个机组取出对应系数
                    for w = 1:W
                        Zon = [Zon,0,zeros(1,bN),Zmid(1,2+(tt-1)*W*bpN+(w-1)*bpN:1+(tt-1)*W*bpN+w*bpN)];
                    end
                end
                cons = [cons, Wc'* Lamda_on == Zon'];
                cons = [cons, Uc'* Lamda_on + Lamda2_on== 0];
                cons = [cons, hc'* Lamda_on >= 0];
                cons = [cons, Lamda2_on >= 0];
            end
            %关机限制
            omiga = max(0,t-Toff(index))+1;
            if (t >= 1+L0(index)&&(omiga<=t)) % 等价于(t >= 1+L0(g))&&(t <= T)，因为循环内t一定小于等于T所以后半部分限制不要了。
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
                for tt = 1:t %按当前时段用到的前t个时段W个机组取出对应系数
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
%         % 线路功率约束
        for i = 1:size(all_branch.I,1) %从第一条支路开始循环遍历所有支路
            left = all_branch.I(i); %支路起点和终点即可得到B电纳
            right = all_branch.J(i);
            abs_x4branch = abs(1/B(left,right));  % 当前支路的阻抗的绝对值 |x_ij|
            % 不等式右端项
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
            % 不等式左端项
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
% end of 不等式约束统一构造

% 等式约束
% 开关机状态等式约束
for g = 1:G
    for t = 1:T
        if t == 1
            Zmot_0 = [0,0,Xs(g,t)] - [0,0,Xd(g,t)] - [0,0,Xu(g,t)];
            Zmot = [0,zeros(1,bN),Xs(g,T+1:T+W*bpN)] - [0,zeros(1,bN),Xd(g,T+1:T+W*bpN)] - [0,zeros(1,bN),Xu(g,T+1:T+W*bpN)]; %无整型变量，且t=1时，uit-1为0；
        else
            Zmot_0 = [0,0,Xs(g,t)-Xd(g,t)-Xu(g,t)+Xu(g,t-1)];
            Bmot = Xs(g,T+1+Sumt(t-1)*W*bpN:T+Sumt(t)*W*bpN) - Xd(g,T+1+Sumt(t-1)*W*bpN:T+Sumt(t)*W*bpN) - Xu(g,T+1+Sumt(t-1)*W*bpN:T+Sumt(t)*W*bpN) + [Xu(g,T+1+Sumt(t-2)*W*bpN:T+Sumt(t-1)*W*bpN),zeros(1,W*bpN)]; %无整型变量，且t!=1时，uit-1要加上；
            Zmot = [];
            for tt = 1:t
                for w = 1:W
                    Zmot = [Zmot,0,zeros(1,bN),Bmot(1,1+(tt-1)*W*bpN+(w-1)*bpN:(tt-1)*W*bpN+w*bpN)];
                end
            end
        end
        cons = [cons, Zmot_0*ones(3,1) == 0]; %等式中无常数项，因此Z0无需考虑
        cons = [cons, Zmot == 0]; %等式中无常数项，因此Z0无需考虑
    end
end
% end of 开关机状态等式约束

% 直流潮流约束等式形式
for n = 1:N %遍历所有节点
    for t = 1:T
        Zac = zeros(1,1+t*W*bN); %连续型变量决策规则
        dac = 0; %常数项
        Zac = Zac + [Yl(n,t),Yl(n,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN)];
        if ismember(n,Gbus) %判断是否是火力机组节点
            index = find(Gbus==n);
            Zac = Zac + [YG(index,t),YG(index,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN)];
        end
        if ismember(n,Wbus) %判断是否是风电机组节点
            index = find(Wbus==n);
            Zac = Zac + [YPw(index,t),YPw(index,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN)];
        end
        for j = 1:N %加上节点潮流
            if B(n,j) ~=0
                Zac = Zac - B(n,j) * [Ytheta(j,t),Ytheta(j,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN)]; 
            end
        end
        if ismember(n,Dbus) %判断是否是负载节点
            index = find(Dbus==n);
            dac = dac + SCUC_data.busLoad.node_P(t,index);
        end
        cons = [cons, Zac(1,1) == dac]; 
        cons = [cons, Zac(2:t*W*bN)' == 0]; 
    end
end
% end of 直流潮流约束等式形式
% end of 等式约束

% 构造目标函数约束
Cobj = zeros(1,T*W); %ksi系数矩阵
Aobj = zeros(1,T*W*bN); %ksi_hat系数矩阵
Bobj = zeros(1,T*W*bpN); %ksi_tine系数矩阵
Cobj0 = zeros(1,T); %ksi0系数矩阵
Aobj0 = zeros(1,T); %ksi_hat0系数矩阵
Bobj0 = zeros(1,T); %ksi_tine0系数矩阵
Uobj_l = []; %目标函数凸包约束等式约束顶点矩阵
Uobj_c = []; %目标函数凸包约束等式约束lamda求和为1系数矩阵
for t = 1:T
    for n = 1:N
        Aobj0 = Aobj0+ 100*[zeros(1,t-1),Yl(n,t),zeros(1,T-t)];
        Aobj = Aobj + 100*[Yl(n,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN),zeros(1,(T-t)*W*bN)]; 
        if ismember(n,Gbus) %判断是否是火力机组节点
            index = find(Gbus==n);
            Aobj0 = Aobj0+ [zeros(1,t-1),Yz(index,t),zeros(1,T-t)] + [zeros(1,t-1),YS(index,t),zeros(1,T-t)];
            Bobj0 = Bobj0 + [zeros(1,t-1),Xs(index,t),zeros(1,T-t)];
            Aobj = Aobj + [Yz(index,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN) + YS(index,T+1+Sumt(t-1)*W*bN:T+Sumt(t)*W*bN), zeros(1,(T-t)*W*bN)]; %将所有zit和Sit决策系数对应加起来
            Bobj = Bobj + [Chot(index) * Xs(index,1+Sumt(t-1)*W*bpN:Sumt(t)*W*bpN), zeros(1,(T-t)*W*bpN)]; %将所有sit决策系数对应加起来
        end
    end 
    Ulamda = blkdiag(verts(:,1+(t-1)*2*vN:vN+(t-1)*2*vN),verts(:,1+vN+(t-1)*2*vN:2*vN*t)); %顶点凸组合表示ksi约束的系数矩阵
    Uconstant = []; %顶点系数lamda组合为1约束的系数矩阵
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
hobj = [zeros(1,T*W*kl+3*T),ones(1,T*W+3*T)]'; %凸包常数项
v = sdpvar(1,M); %构造第一次lagrangian等式约束乘子v
beta = sdpvar(1); %构造第一次lagrangian不等式约束乘子beta
y1 = sdpvar((T*W*(kl+1))+2*3*T,M);%构造第二次等式lagrangian乘子y1，与h同维，每个样本都有一个这样的乘子向量
constraints = [];
constraints = [constraints, beta >= 0]; %不等式乘子要大于等于0
for m = 1:M %开始构造目标函数转换带出的约束
    y_j = y1(:,m); 
    mid_obj = Zobj' - W_obj'*y_j;
    constraints = [constraints, (mid_obj' * ksi_wan{m} + hobj' * y_j) <= M*v(1,m)];
    constraints = [constraints, abs(mid_obj) <= beta];
    constraints = [constraints, U' * y_j >= 0];
end
% end of 构造目标函数约束

constraints = [constraints, cons]; %添加物理约束

objective = sum(v) + eps*beta; %设置目标函数

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

    real_PG = zeros(G,T); %实际火力发电量
    real_z = zeros(G,T); %实际火力发电费用
    real_S = zeros(G,T); %实际火力机组开机费用
    real_l = zeros(N,T); %实际负载切割量
    real_PW = zeros(W,T); %实际风力发电量
    real_theta = zeros(N,T); %实际节点相角
    real_u = zeros(G,T); %实际机组状态
    real_s = zeros(G,T); %实际开机动作
    real_d = zeros(G,T); %实际关机动作
   
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
    for n = 1:N %遍历所有节点
        for t = 1:T
            real_yyy(n,t) = real_yyy(n,t) + real_l(n,t);
            if ismember(n,Gbus) %判断是否是火力机组节点
                index = find(Gbus==n);
                real_yyy(n,t) = real_yyy(n,t) + real_PG(index,t);
            end
            if ismember(n,Wbus) %判断是否是风电机组节点
                index = find(Wbus==n);
                real_yyy(n,t) = real_yyy(n,t) + real_PW(index,t);
            end
            for j = 1:N %加上节点潮流
                if B(n,j) ~=0
                    real_yyy(n,t) = real_yyy(n,t) - B(n,j)*real_theta(j,t);
                end
            end
            if ismember(n,Dbus) %判断是否是负载节点
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