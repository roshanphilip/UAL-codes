%=======================================================================%
% Ref: R.P. Saji, P. Pantidis, M. Mobasher, A new unified arc-length method for damage mechanics problems, manuscript submitted for publication
% Author 1: Roshan Philip Saji
% email: rs7625@nyu.edu
% Author 2: Dr. Panos Pantidis
% email: pp2624@nyu.edu
% Author 3: Prof. Mostafa Mobasher 
% email: mostafa.mobasher@nyu.edu
% Computational Solid Mechanics Lab, New York University Abu Dhabi
% Date: 24-Aug-2023 
%=======================================================================%

clc
clear 
close all
format compact

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%============================ USER INPUTS ================================%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Enter location of current file path
main_file_path='C:\Users\rs7625\Desktop\Unified 2D Codes';
% Enter location of results storage folder 
save_path='C:\Users\rs7625\Desktop\Unified 2D Codes\Saved Files'; 

cd (main_file_path);

% Enter your model name
model_name = "SNS";

% Enter convergence tolerance
tolerance=1e-8;

% For UAL routine:
% Enter the upper limit of ArcLength 
ArcLength_upper_limit = 1e-2;
% Enter the lower limit of ArcLength 
ArcLength_lower_limit = 1e-24;

% I] Local Damage
% a) Choose the solution scheme to be used:
     Scheme_ID = 1;            % 1 - Partitioned Consistent (PC) scheme   2 - Partitioned Non-Consistent (PNC) scheme 
% b) Choose a strain tolerance (ST) value 
     ST = 1e-6;                % Recommended range: 1e-5 to 1e-7
     
% =========================================================================
% ============================== IMPORTANT ===============================%
% =========================================================================
% ------------------------ Parameters file guide --------------------------
% -------------------------------------------------------------------------
% Choose solution scheme 
% SolverID:        1 - Local                (No._DOF_per_node: 2)
%                  2 - Nonlocal Gradient    (No._DOF_per_node: 3)   

% NOTE: Ensure that the "No._DOF_per_node" is updated in the parameter file
% when switching from Local to Nonlocal Gradient model

% TangentID:       1 - Analytical 
%                  2 - Numerical
%      
% RoutineID:       1 - Unified Arclength 
%                  2 - Newton Raphson
%            
% Constraint type: 1 - Cylindrical 
%                  2 - Spherical
% -------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ======================= INCLUDE GLOBAL VARIABLES ========================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

func_include_flags;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ============================ LOAD INPUT FILE ============================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
infile     = fopen(model_name + "_parameters.txt",'r');
[SolverID,TangentID,RoutineID,ncoord,ndof,lc, ... 
            increment,inc_success_counter,min_iter,max_iter,max_accept_iter, ... 
            loadfactor,dlfactor,dlfactor_incr_threshold,increment_plot_threshold,loadfactor_plot_threshold, ...
            flaglf,countflaglf,incrflag,flagplot, ...
            ndomains,nprops,materialprops,alpha_val,beta_val,e_delta,dmax, ...
            nnodes,coords,nelem,maxnodes,connect,nelnodes,elident_vec,nfix,fixnodes,ArcLength_0,Constraint_type,delta_m_bar_0,Applied_Force_Load] = func_read_input_file(infile);  

fclose(infile);

if RoutineID==1
    strain_tolerance = ST;   
else
    strain_tolerance = [];
end

% -------------------------------------------------------------------------
% Find the number of elements attached to each node
num_elem_at_node = zeros(nnodes,1);
for i = 1:size(connect,1)
    for j = 1:size(connect,2)
        num_elem_at_node(connect(i,j)) = num_elem_at_node(connect(i,j)) + 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ================== SPECIFY VALUES FOR GLOBAL VARIABLES ==================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% Non-local gradient parameter (g = lc^2/2)
g = lc^2/2;

% -------------------------------------------------------------------------
% Element length
lelem = coords(1,2) - coords(1,1);

if SolverID == 1
    solver_string = "Local";
elseif SolverID == 2
    solver_string = "Nonlocal_gradient";
elseif SolverID == 3
    solver_string = "Nonlocal_integral";
else
    disp("Check your SolverID - FEM Main Script")
end

if TangentID == 1
    tangent_string = "Analytical";
elseif TangentID == 2
    tangent_string = "Numerical";
else
    disp("Check your TangentID - FEM_Main_Script")
end

if RoutineID == 1
    routine_string = "UAL";
elseif RoutineID == 2
    routine_string = "NR";
elseif RoutineID == 3
    routine_string = "FAL";
else
    disp("Check your RoutineID - FEM_Main_Script")
end

% -------------------------------------------------------------------------
% Calculate matrix with neighbouring points and their weights for the nonlocal integral method
if SolverID == 1 || SolverID == 2
    n_hood = [];
    weights = [];
elseif SolverID == 3
    [n_hood, weights] = func_nhood_gausspts(lc);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====================== SETTING UP MATRICES/VECTORS ======================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Delastic     = func_Delastic;           % Calculate the elastic constitutive matrix D (3x3 for 2D problems)
fixnodes     = (sortrows(fixnodes'))';  % Ensure that the fixnodes matrix lists the prescribed displacements in ascending order of nodes and dofs (as in the global system)
plot_storage = zeros(1,2);              % Saves the converged displacement and reaction forces for all converged increments 

% UAL input initialization
conv_delta_dofs               = zeros(ndof*nnodes,1);                       % Contains converged values of delta_dofs
conv_delta_f_rct_disp_ebc     = zeros(nfix,1);                              % Contain converged values of delta_f_rct_disp_ebc
nfree                         = (ndof*nnodes)-nfix;                         % Contains number of free nodes
Applied_Displacement_Load     = (fixnodes(3,:))';                           % Contains loads applied at the prescribed nodes
delta_m_bar                   = delta_m_bar_0;                              % Initial value of Arc Length load factor
history_var_mat               = zeros(nelem,4);                             % Contains the values of the history variable (damage or kappa) at each Gauss point
[history_var_mat_conv, strain_var_mat_conv]   = deal(zeros(nelem,4));       % Contains the values of the strain variable (damage or kappa) at each Gauss point
[u_bar,u_bar_conv]            = deal(zeros(size(Applied_Displacement_Load,1),1)); % Prescribed displacement vector 
Res_F_F                       = zeros(nfree,1);                                   % Initialize the Residual (Free dofs) vector
Res_F_E_rct                   = zeros(nfix,1);                                    % Initialize the Residual (Prescribed dofs) vector
[convergance_flag,delta_e_nl_conv,e_nl_conv,g_constraint,ArcLength,Load_percentage_last_saved,opt_counter,delta_m_bar_conv,delta_u_bar_conv,delta_u_f_conv,delta_f_rct_disp_ebc_conv,delta_u_bar,delta_u_f,delta_f_rct_disp_ebc,m_bar,u_f,f_rct_disp_ebc,m_bar_conv,u_f_conv,f_rct_disp_ebc_conv]=deal(zeros);

% NR input initialization
[dofs_stored,dofs]          = deal(zeros(ndof*nnodes,1)); % Contains the values of displacements at each node
local_strain_mat_stored     = zeros(nelem,4);       % Contains the values of local equivalent strain at each Gauss point
nonlocal_strain_mat_stored  = zeros(nelem,4);       % Contains the values of nonlocal equivalent strain at each Gauss point
history_var_mat_stored      = zeros(nelem,4);       % Contains the values of the history variable (damage or kappa) at each Gauss point
stress_s1_mat_stored        = zeros(nelem,4);       % Contains the values of 1st principal stress at each Gauss point
xcoord_GP_mat               = zeros(nelem,4);       % Contains the x-coordinates of all Gauss points     
ycoord_GP_mat               = zeros(nelem,4);       % Contains the y-coordinates of all Gauss points
[sum_internal_force,sum_external_force,f_reaction_sum_x,f_reaction_sum_y] =deal(0);

% FAL input initialization
if RoutineID == 3
    % Identify the position of the free and essential nodes
    [~, ~,~,ID_dofs_list_at_ebc,ID_dofs_list_at_nbc,ID_q_dofs,n_essential,n_free]         = func_partitiond_FAL(fixnodes, dofs_stored,nnodes,Applied_Force_Load,[],[],[],1);
    del_u_bar                                                                             = zeros(size(ID_dofs_list_at_ebc,1),1);
    [lambda,lambda_last_converged,delta_lambda_lastconverged,delta_lambda]                = deal(0);
    [u_last_converged,delta_u_last_converged,delta_u]                                     = deal(zeros((size(ID_dofs_list_at_ebc,1)+size(ID_dofs_list_at_nbc,1)),1));
    del_u_f                                                                               = zeros(size(ID_dofs_list_at_nbc,1),1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ======================= GAUSS POINTS COORDINATES ========================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the 4x4 shape function matrix (constant for all our elements)
xilist  = func_integrationpoints;  % Positions of integration points in parent element
for integ_point = 1:4
    xi = xilist(:,integ_point);                 % Size: 2x1
    N(:,integ_point) = func_shapefunctions(xi); % Size: 4x1
    dNdx = func_shapefunctionderivs(xi);

end

% Find and store the x and y coordinates of all Gauss points
for lmn = 1:nelem

lmncoord = zeros(ncoord,maxnodes);
    % Extract nodal coordinates for the current element
    for a = 1:nelnodes(lmn)
        for i = 1:ncoord
            lmncoord(i,a) = coords(i,connect(a,lmn));
        end        
    end
    
    % Transpose lmncoord (from 2x4 to 4x2)
    lmncoord_transp = lmncoord';
    % Calculate the x and y coordinates of the Gauss points
    xcoord_GP_mat(lmn,:) = N' * lmncoord_transp(:,1);
    ycoord_GP_mat(lmn,:) = N' * lmncoord_transp(:,2);

end

% % Identify the position IDs of different nodes based on different
% % classifications
[ID_dofs_list_u_p,ID_dofs_list_u_f,ID_dofs_list_nl_strain,ID_dofs_list_disp,ID_dofs_list_x,ID_dofs_list_y,ID_free_nodes_e,ID_prescribed_nodes_e] = func_find_dof_IDs(fixnodes,ndof,dofs_stored,SolverID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ============================ UAL ANALYSIS ===============================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initially set the residual above the tolerance to start the analysis
residual_norm=2*(tolerance);

if RoutineID == 1
    
    % Calculate the reference dof to measure the displacement applied at each increment   
    [~, ref] = max(Applied_Displacement_Load);

    while abs(u_bar_conv(ref,1)) <= abs(Applied_Displacement_Load(ref,1)) || residual_norm(end)>tolerance
        IsProj = 0;
        % Perform the Unified Arclength analysis
        [convergance_flag,delta_e_nl_conv,e_nl_conv,delta_m_bar,increment,J, dofs, residual_norm, g_constraint,Res_F_E_rct,Res_F_F,last_iteration, history_var_mat_conv, strain_var_mat_conv, ~, ~,Reactions_x,Reactions_y,delta_m_bar_conv,delta_u_bar_conv,delta_u_f_conv,delta_f_rct_disp_ebc_conv,m_bar_conv,u_bar_conv,u_f_conv,f_rct_disp_ebc_conv,ArcLength,Load_percentage_last_saved] = func_UnifiedArcLength_DispControl(model_name,g_constraint,Res_F_E_rct,Res_F_F,Constraint_type,strain_var_mat_conv,tolerance,Load_percentage_last_saved,main_file_path,save_path,ArcLength_0,Applied_Displacement_Load,delta_m_bar_conv,delta_u_bar_conv,delta_u_f_conv,delta_f_rct_disp_ebc_conv,m_bar_conv,u_bar_conv,u_f_conv,f_rct_disp_ebc_conv,ArcLength,delta_m_bar_0,dofs,fixnodes,increment,Delastic,history_var_mat_conv,num_elem_at_node,n_hood,weights,Scheme_ID,strain_tolerance,IsProj,RoutineID,delta_m_bar,ID_dofs_list_u_p,ID_dofs_list_u_f,ID_dofs_list_nl_strain,ID_dofs_list_disp,ID_free_nodes_e,ID_prescribed_nodes_e,delta_e_nl_conv,e_nl_conv,convergance_flag);

        % -----------------------------------------------------------------
        % ------------ SAVING AND ADAPTIVE ARCLENGTH ROUTINE --------------
        % -----------------------------------------------------------------

        % Save convergred increments
        if residual_norm(1,end)<tolerance
            % Assign Force Displacement plot values
            plot_storage(increment+1,1)=abs(u_bar_conv(ref,1));
            plot_storage(increment+1,2)=Reactions_y;                

            % Update increment
            increment=increment+1;

            % Save converged variables for the current iteration
            format shortg
            Load_percentage=round(((u_bar_conv(ref,1)/Applied_Displacement_Load(ref,1))*100),4);
            file_name=sprintf('Increment = %d, Load percentage = %.3f .mat',increment-1,Load_percentage);
            Load_percentage_last_saved=Load_percentage;
            cd (save_path)
            save (file_name)
            cd (main_file_path)
            format compact

            % Update ArcLength
            if last_iteration<5
                ArcLength=min(ArcLength_upper_limit,(10^(log10(ArcLength)+0.2)));           
            elseif last_iteration>12
                ArcLength=max(ArcLength_lower_limit,(10^(log10(ArcLength)-0.2)));
            end
        else
            % Update ArcLength
            ArcLength=max(ArcLength_lower_limit,(10^(log10(ArcLength)-0.2)));
        end

        % Exit if the ArcLength reaches the lower limit
        if ArcLength<=ArcLength_lower_limit
            exit
        end
    end


    % -----------------------------------------------------------------
    % ----------------- PLOTTING PROPERTIES ROUTINE -------------------
    % -----------------------------------------------------------------
    IsProj = 1;
    [~,~,~,~,~,~, ~, ~, ~,~,~,~, ~, ~, gausspoints_prop_mat, nodes_prop_mat,~,~,~,~,~,~,~,~,~,~,~,~] = func_UnifiedArcLength_DispControl(model_name,g_constraint,Res_F_E_rct,Res_F_F,Constraint_type,strain_var_mat_conv,tolerance,Load_percentage_last_saved,main_file_path,save_path,ArcLength_0,Applied_Displacement_Load,delta_m_bar_conv,delta_u_bar_conv,delta_u_f_conv,delta_f_rct_disp_ebc_conv,m_bar_conv,u_bar_conv,u_f_conv,f_rct_disp_ebc_conv,ArcLength,delta_m_bar_0,dofs,fixnodes,increment,Delastic,history_var_mat_conv,num_elem_at_node,n_hood,weights,Scheme_ID,strain_tolerance,IsProj,RoutineID,delta_m_bar,ID_dofs_list_u_p,ID_dofs_list_u_f,ID_dofs_list_nl_strain,ID_dofs_list_disp,ID_free_nodes_e,ID_prescribed_nodes_e,delta_e_nl_conv,e_nl_conv,convergance_flag);
    max_damage = max(gausspoints_prop_mat(:,1));
    func_plotmesh2(coords,connect,nelem,nodes_prop_mat(:,1),nodes_prop_mat(:,8),nodes_prop_mat(:,9),nelnodes,'interp',model_name,inc_success_counter,last_iteration,increment,loadfactor,solver_string,tangent_string,SolverID);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ======================= NEWTON - RAPHSON ANALYSIS =======================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if RoutineID == 2

    % Calculate the reference dof to measure the displacement applied at each increment    
    [ref, ref_node_ID, ref_node_pos] = finc_identify_ref(fixnodes, ndof, model_name);  

    while loadfactor <= 1 && countflaglf == 0 && incrflag < 100

        % Condition to terminate the while loop: the total load is applied
        if flaglf == true; countflaglf = 1; end

        % Condition to terminate the while loop: dlfactor is too small
        if dlfactor < 10^-30; break; end

        IsProj = 0;

        % Perform the Newton Raphson analysis
        [J, dofs, Res_F_F_norm, Res_u_norm, last_iteration, history_var_mat, ~, ~, Reactions_x, Reactions_y] = func_NewtonRaphson_DispControl(model_name,RoutineID,dofs_stored,[fixnodes(1:2,:); fixnodes(3,:)*loadfactor],increment,Delastic,history_var_mat_stored,num_elem_at_node,n_hood,weights,tolerance,IsProj);

        % ---------------------------------------------------------------------
        % ------------------------ ADAPT LOAD ROUTINE -------------------------
        % ---------------------------------------------------------------------

        % Adapt load incrementation w.r.t. the number of iterations needed for convergence
        if last_iteration <= min_iter
            % Fast Convergence - dlfactor increases
            [dofs_stored,loadfactor_stored,history_var_mat_stored,incrflag,flagplot,loadfactor,increment, ...
                inc_success_counter,flaglf,dlfactor] = func_adapt_quick_convergence(dofs,loadfactor,history_var_mat,last_iteration,dlfactor,increment,inc_success_counter);

        elseif (min_iter < last_iteration) && (last_iteration <= max_iter)
            % Moderate Convergence - dlfactor remains the same
            [dofs_stored,loadfactor_stored,history_var_mat_stored,incrflag,flagplot,loadfactor,increment, ...
                inc_success_counter,flaglf] = func_adapt_moderate_convergence(dofs,loadfactor,history_var_mat,dlfactor,increment,inc_success_counter);

        elseif (max_iter < last_iteration) && (last_iteration < max_accept_iter)
            % Slow Convergence - dlfactor decreases
            [dofs_stored,loadfactor_stored,history_var_mat_stored,incrflag,flagplot,loadfactor,increment, ...
                inc_success_counter,flaglf,dlfactor] = func_adapt_slow_convergence(dofs,loadfactor,history_var_mat,last_iteration,dlfactor,increment,inc_success_counter);

        else
            % No Convergence - discard the last step and repeat with smaller load value
            [flagplot,flaglf,countflaglf,incrflag,loadfactor,dlfactor] = func_adapt_no_convergence(incrflag,loadfactor_stored,dlfactor);

        end

        % -----------------------------------------------------------------
        % ----------------------- SAVING ROUTINE --------------------------
        % -----------------------------------------------------------------

        % Save convergred increments
        if Res_u_norm(1,end)<tolerance
            % Assign Force Displacement plot values
            plot_storage(increment,1)=abs(dofs(ref,1));
            plot_storage(increment,2)=Reactions_y;                

            % Save converged variables for the current iteration
            format shortg
            Load_percentage=round(((dofs(ref,1)/fixnodes(3,ref_node_pos))*100),4);
            file_name=sprintf('Increment = %d, Load percentage = %.3f .mat',increment-1,Load_percentage);
            Load_percentage_last_saved=Load_percentage;
            cd (save_path)
            save (file_name)
            cd (main_file_path)
            format compact
        end
       
        % -----------------------------------------------------------------
        % ----------------- PLOTTING PROPERTIES ROUTINE -------------------
        % -----------------------------------------------------------------

        if (last_iteration ~= max_accept_iter) && (loadfactor > loadfactor_plot_threshold) % (rem(increment,increment_plot_threshold) == 0)
            IsProj = 1;
            [~, ~, ~, ~, ~, ~, gausspoints_prop_mat, nodes_prop_mat, ~, ~] = func_NewtonRaphson_DispControl(model_name, RoutineID, dofs, [fixnodes(1:2,:); fixnodes(3,:)*loadfactor], increment, Delastic, history_var_mat_stored, num_elem_at_node, n_hood, weights, tolerance, IsProj);
            max_damage = max(gausspoints_prop_mat(:,1));   
            func_plotmesh2(coords,connect,nelem,nodes_prop_mat(:,1),nodes_prop_mat(:,8),nodes_prop_mat(:,9),nelnodes,'interp',model_name,inc_success_counter,last_iteration,increment,loadfactor,solver_string,tangent_string,SolverID);                    
        end

    end
end