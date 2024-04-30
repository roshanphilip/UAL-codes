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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%============================== USER INPUT ================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Enter location of current file path
file_path = 'C:\Users\rs7625\Desktop\UAL-codes-main\Unified 1D Codes';
% Enter location of results storage folder
save_path = 'C:\Users\rs7625\Desktop\UAL-codes-main\Unified 1D Codes\Saved Files';

% Enter results file name
file_name = sprintf('Example.mat');

% Go to current file path
cd (file_path)

% Plotting triggers
plot_flag=1; % 0 - Don't plot results; 1 - Plot results
save_flag=1; % 0 - Don't save results; 1 - Save results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%============================ RUN INPUT FILE ==============================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input parameters
[delta_dofs_con, ST, ArcLength_0_FAL, max_ArcLength_FAL, min_ArcLength_FAL, max_displacement_FAL, lambda_con, delta_lambda_con, n_FreeNodes, n_EssentialNodes, reset_counter, u_f_con, delta_f_rct_con, f_rct_con, e_nl_con, delta_u_p_con, max_ArcLength,min_ArcLength,f_reaction_essential_con,m_bar_con,delta_u_bar_con, delta_u_f_con, delta_f_reaction_essential_con, delta_m_bar_con, delta_e_nl_con, ArcLength, n_dofs,c, u_p_con, storage_damage_loadstep_con, Applied_Force_Load, plot_switch, max_delta_lf, min_delta_lf, lf, SolverID,l_Total,n_NodesPerElement,weight_integrationpoint,a,b,n_TotalElements,constraint_type,Applied_Displacement_Load,tolerance,k_max,delta_lf,delta_m_bar_0,DamageThreshholdStrain,n_TotalNodes,n_IntegrationPoints,ConstitutiveMatrix,YoungsModulus,Globalcoordinates_EachNode,connect_LocaltoGlobal_Nodenumbers,n_Global,storage_damage_loadstep,plot_storage,dofs_con,RoutineID,f_reaction] = Inputs_1D();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==========================INITIALIZE VARIABLES============================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Identify dof IDs
[ID_prescribed_u,ID_free_u,ID_u,ID_nl_strain,ID_free, q] = finc_ID_dofs(n_dofs, n_TotalNodes, Applied_Force_Load, dofs_con, SolverID, RoutineID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%============================= UAL ANALYSIS ===============================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve the problem based on the chosen aprroximation method
if RoutineID == 1
    % Start measuring computational time
    tic

    % Initialize increment counter
    n=1;

    while max(u_p_con)<Applied_Displacement_Load(2,1)
        % UAL solver
        [delta_m_bar_con,delta_u_p_con,delta_u_f_con,delta_e_nl_con, delta_f_rct_con,m_bar_con, u_p_con, u_f_con, e_nl_con, f_rct_con, n, k, storage_damage_loadstep, storage_strain_loadstep, storage_nonlocal_strain_loadstep, residual_norm, ArcLength] = func_UAL(ST, delta_m_bar_0, constraint_type, n_dofs, c, delta_m_bar_con,delta_u_p_con,delta_u_f_con,delta_e_nl_con, delta_f_rct_con, m_bar_con, u_p_con, u_f_con, e_nl_con, f_rct_con,n_NodesPerElement,weight_integrationpoint,a,b,n_TotalElements,Applied_Displacement_Load,tolerance,k_max,n,ArcLength,DamageThreshholdStrain,n_TotalNodes,n_IntegrationPoints,ConstitutiveMatrix,YoungsModulus,Globalcoordinates_EachNode,connect_LocaltoGlobal_Nodenumbers,n_Global,storage_damage_loadstep_con,ID_prescribed_u,ID_free_u,ID_u,ID_nl_strain,ID_free,SolverID);
        
        % If convergence is reached, save the state and history variables
        if  residual_norm(1,end)<tolerance

            % Save the history variable (damage)
            storage_damage_loadstep_con = storage_damage_loadstep;

            % Save the current damage and strains across the domain
            storage_damage_loadstep_global_matrix(n,:)           =   storage_damage_loadstep;
            storage_strain_loadstep_global_matrix(n,:)           =   storage_strain_loadstep;
            storage_nonlocal_strain_loadstep_global_matrix(n,:)  =   storage_nonlocal_strain_loadstep;

            % Save the F-displacement parameters
            plot_storage(n+1,1)                                  =   max(u_p_con);
            plot_storage(n+1,2)                                  =   max(f_rct_con);

            % Update increment counter
            n = n+1;

            % Update delta_l_f
            if k<5
                ArcLength  = min(max_ArcLength,(10^(log10(ArcLength)+0.2)));
            elseif k>12
                ArcLength  = max(min_ArcLength,(10^(log10(ArcLength)-0.2)));
            end
            % If convergence is not reached, reduce the step size and try again
        else
                ArcLength  = max(min_ArcLength,(10^(log10(ArcLength)-0.2)));

            if (ArcLength <= min_ArcLength)
                disp("REDUCE max_ArcLength")
                disp("===================")
                plot_switch = 0;
                break
            end

        end
    end

    if plot_switch == 1
        % Save computational time taken
        runtime = toc
        if save_flag == 1
            save (file_name)
        end
        % Plotting functions
        func_plot(plot_storage(:,1),plot_storage(:,2),storage_damage_loadstep_global_matrix,storage_strain_loadstep_global_matrix,storage_nonlocal_strain_loadstep_global_matrix,n_TotalElements,l_Total,plot_flag,save_flag,file_path,save_path,SolverID,file_name);
    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%============================= FAL ANALYSIS ===============================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve the problem based on the chosen aprroximation method
if RoutineID == 2
    % Start measuring computational time
    tic

    % Initialize increment counter
    n=1;

    while max(u_f_con)<max_displacement_FAL  % The max. displacement is used only as an exit criteria in FAL
        % FAL solver
        [delta_dofs_con, delta_lambda_con,delta_u_p_con,delta_u_f_con,delta_e_nl_con, delta_f_rct_con,lambda_con, u_p_con, u_f_con, e_nl_con, f_reaction, n, k, storage_damage_loadstep, storage_strain_loadstep, storage_nonlocal_strain_loadstep, residual_norm, ArcLength] = func_FAL(ST, q, n_FreeNodes, n_EssentialNodes, constraint_type, n_dofs, c, delta_lambda_con,delta_u_p_con,delta_u_f_con,delta_e_nl_con, delta_f_rct_con, lambda_con, u_p_con, u_f_con, e_nl_con, f_reaction,n_NodesPerElement,weight_integrationpoint,a,b,n_TotalElements,tolerance,k_max,n,ArcLength,DamageThreshholdStrain,n_TotalNodes,n_IntegrationPoints,ConstitutiveMatrix,YoungsModulus,Globalcoordinates_EachNode,connect_LocaltoGlobal_Nodenumbers,n_Global,storage_damage_loadstep_con,ID_prescribed_u,ID_free_u,ID_u,ID_nl_strain,ID_free,SolverID, delta_dofs_con);

        % If convergence is reached, save the state and history variables
        if  residual_norm(1,end)<tolerance
      
            % Save the history variable (damage)
            storage_damage_loadstep_con = storage_damage_loadstep;

            % Save the current damage and strains across the domain
            storage_damage_loadstep_global_matrix(n,:)           =   storage_damage_loadstep;
            storage_strain_loadstep_global_matrix(n,:)           =   storage_strain_loadstep;
            storage_nonlocal_strain_loadstep_global_matrix(n,:)  =   storage_nonlocal_strain_loadstep;

            % Save the F-displacement parameters
            plot_storage(n+1,1)                                  =   max(u_f_con);
            plot_storage(n+1,2)                                  =   max(f_reaction);

            % Update increment counter
            n = n+1;

            % Update delta_l_f
            if k<5
                ArcLength  = min(max_ArcLength_FAL,(10^(log10(ArcLength)+0.2)));
            elseif k>12
                ArcLength  = max(min_ArcLength_FAL,(10^(log10(ArcLength)-0.2)));
            end
            % If convergence is not reached, reduce the step size and try again
        else
                ArcLength  = max(min_ArcLength_FAL,(10^(log10(ArcLength)-0.2)));

            if (ArcLength <= min_ArcLength_FAL)
                disp("REDUCE max_ArcLength")
                disp("===================")
                plot_switch = 0;
                break
            end

        end
    end

    if plot_switch == 1
        % Save computational time taken
        runtime = toc
        if save_flag == 1
            save (file_name)
        end
        % Plotting functions
        func_plot(plot_storage(:,1),plot_storage(:,2),storage_damage_loadstep_global_matrix,storage_strain_loadstep_global_matrix,storage_nonlocal_strain_loadstep_global_matrix,n_TotalElements,l_Total,plot_flag,save_flag,file_path,save_path,SolverID,file_name);
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%======================== NEWTON-RAPHSON ANALYSIS =========================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Solve the problem based on the chosen aprroximation method
if RoutineID == 3
    % Start measuring computational time
    tic

    % Initialize iteration and increment counters
    k = 0;                                        % Iteration number
    n = 1;                                        % Increment number

    while max(u_p_con)<Applied_Displacement_Load(2,1)
        % Newton-Raphson solver
        [dofs_con,f_reaction,n,k,storage_damage_loadstep,storage_strain_loadstep,storage_nonlocal_strain_loadstep,residual_norm,delta_lf] = func_NR(ST, n_dofs,c, dofs_con, lf,n_NodesPerElement,weight_integrationpoint,a,b,n_TotalElements,Applied_Displacement_Load,tolerance,k_max,n,delta_lf,DamageThreshholdStrain,n_TotalNodes,n_IntegrationPoints,f_reaction,ConstitutiveMatrix,YoungsModulus,Globalcoordinates_EachNode,connect_LocaltoGlobal_Nodenumbers,n_Global,storage_damage_loadstep_con,ID_prescribed_u,ID_free_u,ID_u,ID_nl_strain,ID_free,SolverID);

        % If convergence is reached, save the state and history variables
        if  residual_norm(1,end)<tolerance
            
            % Save the history variable (damage)
            storage_damage_loadstep_con = storage_damage_loadstep;

            % Unroll dofs
            u_p_con                     = dofs_con(ID_prescribed_u,1);

            % Save the current damage and strains across the domain
            storage_damage_loadstep_global_matrix(n,:)           =   storage_damage_loadstep;
            storage_strain_loadstep_global_matrix(n,:)           =   storage_strain_loadstep;
            storage_nonlocal_strain_loadstep_global_matrix(n,:)  =   storage_nonlocal_strain_loadstep;

            % Save the F-displacement parameters
            plot_storage(n+1,1)                                  =   max(u_p_con);
            plot_storage(n+1,2)                                  =   max(f_reaction);

            % Update increment counter
            n=n+1;         
            
            % Update delta_l_f    
            if max(u_p_con) > 0 && max(u_p_con) < 0.01    % If the NR routine needs to take smaller steps in a region to better ensure convergence, enter its bounds here
                if k<5
                    delta_lf  = min(max_delta_lf,(10^(log10(delta_lf)+0.2)));
                elseif k>12
                    delta_lf  = max(min_delta_lf,(10^(log10(delta_lf)-0.2)));
                end
            else
                delta_lf = 1e-3;
            end
        % If convergence is not reached, reduce the step size and try again 
        else  
            delta_lf      = max(min_delta_lf,(10^(log10(delta_lf)-0.2)));
            
            % In highly non-linear local damage problems, increasing the
            % delta_lf value once it reaches the lower limit, helps reach
            % convergence
            if (delta_lf <= min_delta_lf) && reset_counter > 0               
                delta_lf      = min(1e-3,10^(-reset_counter));
                reset_counter = reset_counter - 1;
            elseif  reset_counter <= 0  
                disp("Change NR parameters to achieve convergence")
                disp("==============================================")
                plot_switch = 0;
                break
            end

        end
    end

    if plot_switch == 1
        % Save computational time taken
        runtime = toc
        if save_flag == 1
            save (file_name)
        end
        % Plotting functions
        func_plot(plot_storage(:,1),plot_storage(:,2),storage_damage_loadstep_global_matrix,storage_strain_loadstep_global_matrix,storage_nonlocal_strain_loadstep_global_matrix,n_TotalElements,l_Total,plot_flag,save_flag,file_path,save_path,SolverID,file_name);
    end
    

end

