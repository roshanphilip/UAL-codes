%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ======================= DEFINE GLOBAL VARIABLES =========================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FEM MODEL PARAMETERS 
global nprops materialprops ncoord ndof nnodes coords nelem maxnodes connect_nds nelnodes elident_vec nfix fixnodes ndof2
ndof2 = 2;

% NONLOCAL GRADIENT PARAMETER 
% g = lc^2/2, lc = characteristic length
global g

% MAZAR'S DAMAGE MODEL PARAMETERS
global alpha_val beta_val e_delta dmax

% ADAPTIVE LOAD AND PLOTTING PARAMETERS
global min_iter max_iter max_accept_iter dlfactor_incr_threshold increment_plot_threshold loadfactor_plot_threshold

% SOLVER SCHEME
% SolverID: 1 - Local, 2 - Nonlocal Gradient, 3 - Nonlocal Integral
% TangentID: 1 - Analytical, 2 - Numerical
% RoutineID 1 - UAL, 2 - NR 
global SolverID TangentID

% IDs of dofs based on different classifications 
% global ID_dofs_list_u_p ID_dofs_list_u_f ID_dofs_list_nl_strain ID_dofs_list_disp ID_dofs_list_x ID_dofs_list_y
