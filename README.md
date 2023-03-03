# DC-based Security Constraints Formulation: A Perspective of Primal-Dual Interior Point Method

This repository contains the source code of the paper "*DC-based Security Constraints Formulation: A Perspective of Primal-Dual Interior Point Method*"

## Prerequisite

- **MATLAB** to run the code. 
- **YALMIP** to model the optimization problem.
- **Matpower** to get the power system cases.
- **Gurobi** (or other solvers) to solve the optimization problem.

## File list

- `main.m`: The main program for numerical cases. Other functions will be called in this file.
- `OPF_draw_density.m` and `OPF_pure_draw_density.m`: Draw the non-zero elements distribution in the coefficient
matrices in Figure 1.
- `OPF_pure.m` and `SCED.m`: Comparison of size, density and solving time of three formulations for solving the OPF and SCED problem.
- `SCED_del.m` and `SCED_del_v2.m`: Comparison of outer approximations of the SCED problem for both the PTDF and mixed formulations. These two files are 1) randomly removing uncongested branch constraints, and 2) sequentially removing the most uncongested branch constraints, respectively.
- `SCUC.m`: SCUC cases.
