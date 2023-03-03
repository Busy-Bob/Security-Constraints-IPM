function ret = SCUC(mpc, solver)
%%%%
% @para
% - mpc: matpower case
% - solver: string name
% @return 
% - ret.PTDF / ret.theta:       result structure for PTDF method or theta method
% - ret.PTDF.offline_time:      offline aux matrix generation time
% - ret.PTDF.formulation_time:  yalmip formulation time
% - ret.PTDF.solve_time:        solve time (including time of sending data to solver interface)
% - ret.PTDF.yalmip_result:     yalmip result structure
%%%%
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
   MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
   QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

if any(mpc.bus(:, BUS_I) ~= (1:size(mpc.bus, 1))')
    warning("Buses must be numbered consecutively in bus matrix; Start transformation")
    mpc = ext2int(mpc);
end
gen_idx = find(mpc.gen(:, PMAX) > 0); % 
branch_idx = find(mpc.branch(:, BR_STATUS) == 1 & mpc.branch(:, RATE_A) ~= 0);
slack_idx = find(mpc.bus(:, BUS_TYPE) == REF);
gen_bus_idx = mpc.gen(gen_idx, GEN_BUS);
% check case availability (simple)
if sum(1.1 * mpc.bus(:, PD)) > sum(mpc.gen(gen_idx, PMAX))
    error("Total demand is bigger than maximum generation power. Please check the case file.")
end
if size(branch_idx, 1) == 0
    error("Branch RATE_A or BR_STATUS property missed. Please check the case file.")
end
if size(find(mpc.branch(:, BR_STATUS) == 1), 1) ~= size(branch_idx, 1)
    warning("Some branches RATE_A property missed. Set them to infinity.")
end

% parameters
T = 24;
reserve_ratio = 0.1;
UT = 4; 
DT = 4;

% random cost coeff
rng(1); % fix a seed
a = rand(size(gen_idx)) * 0.1;
b = rand(size(gen_idx));
c = rand(size(gen_idx)) * 10;
C_SU = rand(size(gen_idx)) * 50;
C_SD = rand(size(gen_idx)) * 50;
% random power demand coeff per hour
coeff = rand(T, 1) * 0.1 + 0.95; % 0.95~1.05
% random ramping 
RD = mpc.gen(gen_idx, PMAX) .* rand(size(gen_idx)) * 0.2; % 0~0.2
RU = RD;
SD = max(mpc.gen(gen_idx, PMIN), RD); 
SU = SD;

% transform
a = repmat(a, [1 T]);
b = repmat(b, [1 T]);
c = repmat(c, [1 T]);
C_SU = repmat(C_SU, [1 T]);
C_SD =  repmat(C_SD, [1 T]);

P_min = repmat(mpc.gen(gen_idx, PMIN), [1 T]);
P_max = repmat(mpc.gen(gen_idx, PMAX), [1 T]);
branch_max = repmat(mpc.branch(branch_idx, RATE_A), [1 T]);
SD = repmat(SD, [1 T-1]);
SU = repmat(SU, [1 T-1]);
RD = repmat(RD, [1 T-1]);
RU = repmat(RU, [1 T-1]);
P_d = mpc.bus(:, PD) * coeff';
R = sum(P_d, 1) * reserve_ratio;
%% use PTDF
yalmip('clear')
%%%% PTFD generation time start %%%%
tic
% PTDF 
PTDF = makePTDF(mpc.baseMVA, mpc.bus, mpc.branch);
ret.PTDF.offline_time = toc;
%%%% PTFD generation time end %%%%

%%%% formulation time start %%%%
tic
P_g = sdpvar(size(gen_idx, 1), T, 'full');
P_g_bar = sdpvar(size(gen_idx, 1), T, 'full');
r_g = sdpvar(size(gen_idx, 1), T, 'full');
v = binvar(size(gen_idx, 1), T, 'full');
u = binvar(size(gen_idx, 1), T, 'full');
w = binvar(size(gen_idx, 1), T, 'full');
sum_v = [];
sum_w = [];
for k = 1:(T-UT+1)
    sum_v = [sum_v, sum(v(:, k:(k+UT-1)), 2)];
end
for k = 1:(T-DT+1)
    sum_w = [sum_w, sum(w(:, k:(k+DT-1)), 2)];
end
% other constraints
constraints = [];
constraints = [constraints;
    % other constraints
    sum(P_g, 1) == sum(P_d, 1);
    sum(r_g, 1) >= R;
    u(:, 2:T) - u(:, 1:(T-1)) == v(:, 2:T) - w(:, 2:T);
    sum_v <= u(:, UT:T);
    sum_w <= 1 - u(:, UT:T);
    r_g >= 0;
    P_g_bar == r_g + P_g;
    P_min .* u <= P_g <= P_max .* u;
    P_min .* u <= P_g_bar <= P_max .* u;
    P_g_bar(:, 1:T-1) <= P_max(:, 1:T-1) .* (u(:, 1:T-1) - w(:, 2:T)) + SD .* w(:, 2:T);
    P_g_bar(:, 2:T) - P_g(:, 1:T-1) <= RU .* u(:, 1:T-1) + SU .* v(:, 2:T);
    P_g(:, 1:T-1) - P_g(:, 2:T) <= RD .* u(:, 2:T) + SD .* w(:, 2:T);
    % security constraints
    -branch_max <= PTDF(branch_idx, gen_bus_idx) * P_g - PTDF(branch_idx, :) * P_d <= branch_max
    ];

% objective
objective = sum(sum(C_SU .* v + C_SD .* w + a .* P_g .* P_g + b .* P_g + c .* u));
ret.PTDF.formulation_time = toc;
%%%% formulation time end %%%%

%%%% solve time start %%%%
tic
options = sdpsettings('solver', solver, 'verbose', 1);
% options.gurobi.Method = 2;
% options.gurobi.NodeMethod = 2;
ret.PTDF.yalmip_result = optimize(constraints, objective, options);
ret.PTDF.solve_time = toc;
%%%% solve time end %%%%

%% use theta
yalmip('clear')
%%%% aux matrix generation time start %%%%
tic
[Bbus, Bf, Pbusinj, Pfinj] = makeBdc(mpc.baseMVA, mpc.bus, mpc.branch);
C = sparse(gen_bus_idx, 1:size(gen_bus_idx, 1), ones(size(gen_bus_idx)), size(mpc.bus, 1), size(gen_bus_idx, 1)); % bus-generator incidence matrix
ret.theta.offline_time = toc;
%%%% aux matrix generation time end %%%%


%%%% formulation time start %%%%
tic
theta = sdpvar(size(mpc.bus, 1), T, 'full');
P_g = sdpvar(size(gen_idx, 1), T, 'full');
P_g_bar = sdpvar(size(gen_idx, 1), T, 'full');
r_g = sdpvar(size(gen_idx, 1), T, 'full');
v = binvar(size(gen_idx, 1), T, 'full');
u = binvar(size(gen_idx, 1), T, 'full');
w = binvar(size(gen_idx, 1), T, 'full');
sum_v = [];
sum_w = [];
for k = 1:(T-UT+1)
    sum_v = [sum_v, sum(v(:, k:(k+UT-1)), 2)];
end
for k = 1:(T-DT+1)
    sum_w = [sum_w, sum(w(:, k:(k+DT-1)), 2)];
end
% other constraints
constraints = [];
constraints = [constraints;
    % other constraints
    sum(P_g, 1) == sum(P_d, 1);
    sum(r_g, 1) >= R;
    u(:, 2:T) - u(:, 1:(T-1)) == v(:, 2:T) - w(:, 2:T);
    sum_v <= u(:, UT:T);
    sum_w <= 1 - u(:, UT:T);
    r_g >= 0;
    P_g_bar == r_g + P_g;
    P_min .* u <= P_g <= P_max .* u;
    P_min .* u <= P_g_bar <= P_max .* u;
    P_g_bar(:, 1:T-1) <= P_max(:, 1:T-1) .* (u(:, 1:T-1) - w(:, 2:T)) + SD .* w(:, 2:T);
    P_g_bar(:, 2:T) - P_g(:, 1:T-1) <= RU .* u(:, 1:T-1) + SU .* v(:, 2:T);
    P_g(:, 1:T-1) - P_g(:, 2:T) <= RD .* u(:, 2:T) + SD .* w(:, 2:T);
    % security constraints
    Bbus * theta == C * P_g - P_d;
    theta(slack_idx, :) == 0;
    -branch_max <= Bf(branch_idx, :) * theta <= branch_max
    ];

% objective
objective = sum(sum(C_SU .* v + C_SD .* w + a .* P_g .* P_g + b .* P_g + c .* u));
ret.theta.formulation_time = toc;
%%%% formulation time end %%%%

%%%% solve time start %%%%
tic
options = sdpsettings('solver', solver, 'verbose', 1);
% options.gurobi.Method = 2;
% options.gurobi.NodeMethod = 2;
ret.theta.yalmip_result = optimize(constraints, objective, options);
ret.theta.solve_time = toc;
%%%% solve time end %%%%

%% display
disp("1: PTDF offline matrix generation time: " + ret.PTDF.offline_time + "s")
disp("1: PTDF formulation time: " + ret.PTDF.formulation_time + "s")
disp("1: PTDF solve time: " + ret.PTDF.solve_time + "s")
disp("2: theta offline matrix generation time: " + ret.theta.offline_time + "s")
disp("2: theta formulation time: " + ret.theta.formulation_time + "s")
disp("2: theta solve time: " + ret.theta.solve_time + "s")
end