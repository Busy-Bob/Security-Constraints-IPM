clear all; close all; clc;

%% 1. structure
mpc = loadcase(case39); 
OPF_pure_draw_density(mpc, 'gurobi');
OPF_draw_density(mpc, 'gurobi');

%% 2. compare opf
mpc = loadcase(case_ACTIVSg25k); 
options = sdpsettings('solver', 'gurobi', 'verbose', 1);
ret1 = OPF_pure(mpc, options);

mpc = loadcase(case_ACTIVSg10k); 
options = sdpsettings('solver', 'gurobi', 'verbose', 1);
ret1 = OPF_pure(mpc, options);

mpc = loadcase(case6515rte); 
options = sdpsettings('solver', 'gurobi', 'verbose', 1);
ret1 = OPF_pure(mpc, options);

mpc = loadcase(case6470rte); 
options = sdpsettings('solver', 'gurobi', 'verbose', 1);
ret1 = OPF_pure(mpc, options);

mpc = loadcase(case1951rte); 
options = sdpsettings('solver', 'gurobi', 'verbose', 1);
ret1 = OPF_pure(mpc, options);

% without Presolve
mpc = loadcase(case_ACTIVSg25k); 
options = sdpsettings('solver', 'gurobi', 'verbose', 1);
options.gurobi.Presolve = 0;
ret1 = OPF_pure(mpc, options);

mpc = loadcase(case_ACTIVSg10k); 
options = sdpsettings('solver', 'gurobi', 'verbose', 1);
options.gurobi.Presolve = 0;
ret1 = OPF_pure(mpc, options);

mpc = loadcase(case6515rte); 
options = sdpsettings('solver', 'gurobi', 'verbose', 1);
options.gurobi.Presolve = 0;
ret1 = OPF_pure(mpc, options);

mpc = loadcase(case6470rte); 
options = sdpsettings('solver', 'gurobi', 'verbose', 1);
options.gurobi.Presolve = 0;
ret1 = OPF_pure(mpc, options);

mpc = loadcase(case1951rte); 
options = sdpsettings('solver', 'gurobi', 'verbose', 1);
options.gurobi.Presolve = 0;
ret1 = OPF_pure(mpc, options);

% 3. compare SCED
mpc = loadcase(case1354pegase); 
options = sdpsettings('solver', 'gurobi', 'verbose', 1);
ret1 = SCED(mpc, options);

mpc = loadcase(case1951rte); 
options = sdpsettings('solver', 'gurobi', 'verbose', 1);
ret1 = SCED(mpc, options);

mpc = loadcase(case3375wp); 
options = sdpsettings('solver', 'gurobi', 'verbose', 1);
ret1 = SCED(mpc, options);

mpc = loadcase(case6470rte); 
options = sdpsettings('solver', 'gurobi', 'verbose', 1);
ret1 = SCED(mpc, options);

mpc = loadcase(case6515rte); 
options = sdpsettings('solver', 'gurobi', 'verbose', 1);
ret1 = SCED(mpc, options);


%% 4. remove constraints
ratios = [0.8, 0.5, 0.3, 0.1, 0.05, 0.02, 0.01];
mpcs = [case1354pegase, case1951rte, case3375wp, case6470rte, case6515rte];
ptdf_time = zeros(length(mpcs), length(ratios));
theta_time = zeros(length(mpcs), length(ratios));

options = sdpsettings('solver', 'gurobi', 'verbose', 0);
iters = 5; % for average;
for i = 1:length(mpcs)
    for j = 1:length(ratios)
        mpc = mpcs(i);
        ratio = ratios(j);
        disp("-------------------------------------------")
        disp("- mpc: " + size(mpc.bus, 1) )
        disp("- ratio: " + ratio )
        for k = 1:iters
            ret1 = SCED_del(mpc, options, ratio);
        	ptdf_time(i, j) = ptdf_time(i, j) + ret1.PTDF.solve_time / iters;
            theta_time(i, j) = theta_time(i, j) + ret1.theta.solve_time / iters;
        end
    end
end

save('SCED_del_1007.mat', 'mpcs', 'ratios', 'ptdf_time', 'theta_time')

% 
ratios = [0.8, 0.5, 0.3, 0.1, 0.05, 0.02, 0.01];
mpcs = [case1354pegase, case1951rte, case3375wp, case6470rte, case6515rte];
ptdf_time = zeros(length(mpcs), length(ratios));
theta_time = zeros(length(mpcs), length(ratios));

options = sdpsettings('solver', 'gurobi', 'verbose', 0);
iters = 5; % for average;
for i = 1:length(mpcs)
    for j = 1:length(ratios)
        mpc = mpcs(i);
        ratio = ratios(j);
        disp("-------------------------------------------")
        disp("- mpc: " + size(mpc.bus, 1) )
        disp("- ratio: " + ratio )
        for k = 1:iters
            ret1 = SCED_del_v2(mpc, options, ratio);
        	ptdf_time(i, j) = ptdf_time(i, j) + ret1.PTDF.solve_time / iters;
            theta_time(i, j) = theta_time(i, j) + ret1.theta.solve_time / iters;
        end
    end
end
save('SCED_del_1010.mat', 'mpcs', 'ratios', 'ptdf_time', 'theta_time')

