function ret = OPF(mpc, solver)
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
gen_idx = find(mpc.gen(:, GEN_STATUS) > 0);
branch_idx = find(mpc.branch(:, BR_STATUS) == 1 & mpc.branch(:, RATE_A) ~= 0);
slack_idx = find(mpc.bus(:, BUS_TYPE) == REF);
gen_bus_idx = mpc.gen(gen_idx, GEN_BUS);
% check case availability (simple)
if sum(mpc.bus(:, PD)) < sum(mpc.gen(gen_idx, PMIN)) 
    error("Total demand is smaller than minimum generation power. Please check the case file.")
end
if sum(mpc.bus(:, PD)) > sum(mpc.gen(gen_idx, PMAX))
    error("Total demand is bigger than maximum generation power. Please check the case file.")
end
if size(branch_idx, 1) == 0
    error("Branch RATE_A or BR_STATUS property missed. Please check the case file.")
end
if size(find(mpc.branch(:, BR_STATUS) == 1), 1) ~= size(branch_idx, 1)
    warning("Some branches RATE_A property missed. Set them to infinity.")
end

% random cost coeff
rng(1); % fix a seed
a = rand(size(gen_idx)) * 0.1;
b = rand(size(gen_idx));
c = rand(size(gen_idx)) * 10;
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
P_g = sdpvar(size(gen_idx, 1), 1, 'full');
% other constraints
constraints = [];
constraints = [constraints;
    sum(P_g) == sum(mpc.bus(:, PD));
    mpc.gen(gen_idx, PMIN) <= P_g <= mpc.gen(gen_idx, PMAX);
    ];
% security constraints
constraints = [constraints;
    -mpc.branch(branch_idx, RATE_A) <= PTDF(branch_idx, gen_bus_idx) * P_g - PTDF(branch_idx, :) * mpc.bus(:, PD) <= mpc.branch(branch_idx, RATE_A);
    ];
% objective
objective = sum(a .* P_g .* P_g + b .* P_g + c);
ret.PTDF.formulation_time = toc;
%%%% formulation time end %%%%

%%%% solve time start %%%%
tic
options = sdpsettings('solver', solver, 'verbose', 1);
% ret.PTDF.yalmip_result = optimize(constraints, objective, options);
[model,recoverymodel] = export(constraints, objective, options);
figure;
spy(model.A, 'black')
figure;
spy(model.A * model.A' , 'black')

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
P_g = sdpvar(size(gen_idx, 1), 1, 'full');
theta = sdpvar(size(mpc.bus, 1), 1, 'full');
% other constraints
constraints = [];
constraints = [constraints;
    sum(P_g) == sum(mpc.bus(:, PD));
    mpc.gen(gen_idx, PMIN) <= P_g <= mpc.gen(gen_idx, PMAX);
    ];
% security constraints
constraints = [constraints;
    Bbus * theta == C * P_g - mpc.bus(:, PD);
    theta(slack_idx) == 0;
    -mpc.branch(branch_idx, RATE_A) <= Bf(branch_idx, :) * theta <= mpc.branch(branch_idx, RATE_A)
    ];
% 
objective = sum(a .* P_g .* P_g + b .* P_g + c);
ret.theta.formulation_time = toc;
%%%% formulation time end %%%%

%%%% solve time start %%%%
tic
options = sdpsettings('solver', solver, 'verbose', 1);
% ret.theta.yalmip_result = optimize(constraints, objective, options);
[model,recoverymodel] = export(constraints, objective, options);
figure;
spy(model.A, 'black')
figure;
spy(model.A * model.A' , 'black')
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