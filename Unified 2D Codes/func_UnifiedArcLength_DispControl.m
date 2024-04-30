function [convergance_flag,delta_e_nl_conv,e_nl_conv,delta_m_bar,increment,J, dofs_stored, residual_norm, g_constraint,Res_F_E_rct,Res_F_F,iteration, history_var_mat_conv, strain_var_mat_conv, gausspoints_prop_mat, nodes_prop_mat,Reactions_x, Reactions_y,delta_m_bar_conv,delta_u_bar_conv,delta_u_f_conv,delta_f_rct_disp_ebc_conv,m_bar_conv,u_bar_conv,u_f_conv,f_rct_disp_ebc_conv,ArcLength,Load_percentage_last_saved] = func_UnifiedArcLength_DispControl(model_name,g_constraint,Res_F_E_rct,Res_F_F,Constraint_type,strain_var_mat_conv,tolerance,Load_percentage_last_saved,main_file_path,save_path,ArcLength_0,Applied_Displacement_Load,delta_m_bar_conv,delta_u_bar_conv,delta_u_f_conv,delta_f_rct_disp_ebc_conv,m_bar_conv,u_bar_conv,u_f_conv,f_rct_disp_ebc_conv,ArcLength,delta_m_bar_0,dofs_stored,fixnodes_applied,increment,Delastic,history_var_mat_conv,num_elem_at_node,n_hood,weights,Scheme_ID,strain_tolerance,IsProj,RoutineID,delta_m_bar,ID_dofs_list_u_p,ID_dofs_list_u_f,ID_dofs_list_nl_strain,ID_dofs_list_disp,ID_free_nodes_e,ID_prescribed_nodes_e,delta_e_nl_conv,e_nl_conv,convergance_flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ======================= UNIFIED ARCLENGTH METHOD ========================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Include global variables
func_include_flags;

% Recover the last converged values of dofs
dofs=dofs_stored;

% Identify the position of the free and prescribed nodes
[~, u_bar, u_f,~,~,~] = func_partitiond(fixnodes_applied, dofs);

% Create copies of the kappa matrix from the end of the previous increment
history_var_mat_previousinc = history_var_mat_conv;
strain_var_mat_previousinc  = strain_var_mat_conv;

% For local damage
if SolverID == 1 
    [delta_e_nl_conv,delta_e_nl_current,delta_e_nl,e_nl_conv,e_nl,del_e_nl]= deal([]);
end

% For non-local gradient damage
if SolverID == 2
    M_bar_u = size(ID_dofs_list_u_p,1);
    M_u     = size(ID_dofs_list_u_f,1);
    M_bar_e = size(ID_prescribed_nodes_e,1);
    M_e     = size(ID_free_nodes_e,1);
end


% -------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SOLVER ROUTINE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
if IsProj == 0

    % Create empty entries for element/nodal properties
    gausspoints_prop_mat = [];
    nodes_prop_mat       = [];

    % Initiate the iteration counter within each loading step
    iteration = 0;

    % Set residual norm to a high value
    residual_norm=2*tolerance;

    % Unroll dofs vector
    [~, u_bar, u_f,~,ID_dofs_list_at_ebc,ID_dofs_list_at_nbc] = func_partitiond(fixnodes_applied, dofs);

    % Calculate Beta
    [Beta] = func_calc_Beta(Constraint_type,model_name,dofs,Delastic,history_var_mat_previousinc,num_elem_at_node,n_hood,weights,strain_tolerance,strain_var_mat_previousinc,IsProj);

    while residual_norm>tolerance

        % Update iteration counter
        iteration = iteration+1;

        % Evaluate consistent tangent matrix [J] and history variable matrix at the start of the increment
        [~,J, ~, ~, ~, ~, ~, ~,~,~] = func_globalstiffness(model_name,dofs,Delastic,history_var_mat_previousinc,num_elem_at_node,n_hood,weights,strain_tolerance,strain_var_mat_previousinc,IsProj,RoutineID);

        % Partition the consistent tangent matrix
        [J_E, J_EF, J_FE, J_F, ~, ~, ~, ~] = func_partitionK(fixnodes_applied, J);

        if SolverID == 1 % Local damage
            % Calculate corrector coefficients del_u_NR and del_u_g
            del_u_NR=-J_F\Res_F_F;
            del_u_g=-J_F\(J_FE*Applied_Displacement_Load);
        end

        if iteration == 1
            % Calculate predictor values 
            if SolverID == 1  % Local damage
                if increment == 1
                    % Calculate delta predictors %
                    [delta_m_bar,delta_u_bar,delta_u_f,delta_f_rct_disp_ebc,ArcLength_0,ArcLength] = func_delta_predictors_UAL(J_E,J_EF,J_FE,J_F,Applied_Displacement_Load,Beta,del_u_g,del_u_NR,increment,Res_F_E_rct,Res_F_F,g_constraint,delta_u_bar_conv,delta_u_f_conv,delta_f_rct_disp_ebc_conv,delta_m_bar_conv,ArcLength_0,ArcLength);
                else
                    % Calculate delta predictors %
                    [delta_m_bar,delta_u_bar,delta_u_f,delta_f_rct_disp_ebc,ArcLength_0,ArcLength] = func_delta_predictors_UAL(J_E,J_EF,J_FE,J_F,Applied_Displacement_Load,Beta,del_u_g,del_u_NR,increment,Res_F_E_rct,Res_F_F,g_constraint,delta_u_bar_conv,delta_u_f_conv,delta_f_rct_disp_ebc_conv,delta_m_bar_conv,ArcLength_0,ArcLength);
                end
            elseif SolverID == 2 % Non-local gradient damage
                if increment == 1
                     % Paritition Consistent Tangent matrix J
                    [Jpf,Jpp,Jff,Jfp,Jue_p,Jue_f,Jeu_p,Jeu_f,Jee] = func_sort_stiffness_UAL_gradient(J,ID_dofs_list_u_p,ID_dofs_list_u_f,ID_dofs_list_nl_strain,ID_dofs_list_disp,fixnodes_applied);           
                    [delta_m_bar,delta_u_bar,delta_u_f,delta_f_rct_disp_ebc,delta_e_nl,ArcLength] = func_delta_predictors_UAL_gradient(Jff,Jpp,Jpf,Jfp,Jue_p,Jue_f,Jeu_p,Jeu_f,Jee,Applied_Displacement_Load,delta_m_bar_0);
                else
                    [delta_m_bar,delta_u_bar,delta_u_f,delta_f_rct_disp_ebc,delta_e_nl] = func_reset_deltas_UAL_gradient(convergance_flag,delta_m_bar,delta_m_bar_conv,delta_u_bar_conv,delta_u_f_conv,delta_f_rct_disp_ebc_conv,delta_e_nl_conv,Applied_Displacement_Load);
                end
            end
        else
            % Calculate corrector values
            if SolverID == 1 % Local damage
                if Scheme_ID==1
                    [del_f_rct_disp_ebc, del_u_f, del_m_bar,exit_counter] = func_solve_system_equations_consistent_partitioned(J_F,J_E,J_EF,J_FE,Res_F_E_rct,Res_F_F,g_constraint,delta_u_bar,delta_u_f,delta_f_rct_disp_ebc,ID_dofs_list_at_ebc,ID_dofs_list_at_nbc,Beta,Applied_Displacement_Load);
                elseif Scheme_ID==2
                    [del_f_rct_disp_ebc, del_u_f, del_m_bar,exit_counter] = func_solve_system_equations_nonconsistent_partitioned(delta_u_bar_current,delta_u_f_current,delta_f_rct_disp_ebc_current,J_F,J_E,J_EF,J_FE,Res_F_E_rct,Res_F_F,Beta,Applied_Displacement_Load,ArcLength);
                end
            elseif SolverID == 2 % Non-local gradient damage
                [del_f_rct_disp_ebc, del_u_f, del_m_bar,del_e_nl,exit_counter] = func_solve_system_partitioned_UAL_gradient(struct_sub_Matrix_1,Res_F_E_rct,Res_F_F_u_f,g_constraint,Res_F_F_e_star,M_bar_e,M_e);
            end

            % If imaginary roots are present, go back and decrease the ArcLength
            if exit_counter==1
                Restore_arc_length=ArcLength;  
                Restore_residual_norm = residual_norm;
                cd (save_path)
                format shortg
                last_saved_loadstep=sprintf('Increment = %d, Load percentage = %.3f .mat',increment-1,Load_percentage_last_saved);
                load (last_saved_loadstep)
                cd (main_file_path)
                ArcLength=Restore_arc_length; 
                residual_norm = Restore_residual_norm;
                return
            end

            % Update deltas
            delta_m_bar=delta_m_bar+del_m_bar;
            delta_u_bar=delta_m_bar*Applied_Displacement_Load;
            delta_u_f=delta_u_f+del_u_f;
            delta_f_rct_disp_ebc=delta_f_rct_disp_ebc+del_f_rct_disp_ebc;
            
            if SolverID == 2
                delta_e_nl=delta_e_nl+del_e_nl;
            end

        end

        % Save current deltas
        delta_m_bar_current=delta_m_bar;
        delta_u_bar_current=delta_u_bar;
        delta_u_f_current=delta_u_f;
        delta_f_rct_disp_ebc_current=delta_f_rct_disp_ebc;
        
        if SolverID == 2
            delta_e_nl_current=delta_e_nl;
        end

        % Update state variables
        m_bar=m_bar_conv+delta_m_bar_current;
        u_bar=m_bar*Applied_Displacement_Load;
        u_f=u_f_conv+delta_u_f_current;
        f_rct_disp_ebc=f_rct_disp_ebc_conv+delta_f_rct_disp_ebc_current;
       
        if SolverID == 2
            e_nl=delta_e_nl_conv+delta_e_nl_current;
        end

        % Assemble dofs vector
        if SolverID == 1
            dofs(ID_dofs_list_at_ebc) = u_bar;
            dofs(ID_dofs_list_at_nbc) = u_f;
        elseif SolverID == 2
            dofs(ID_dofs_list_u_p) = u_bar;
            dofs(ID_dofs_list_u_f) = u_f;
            dofs(ID_dofs_list_nl_strain)=e_nl;
        end


        % Evaluate consistent tangent matrix [J] and history variable matrix at the current iteration
        [~,J, ~, Res_F, ~,~,history_var_mat,~,~,strain_var_mat] = func_globalstiffness(model_name,dofs,Delastic,history_var_mat_previousinc,num_elem_at_node,n_hood,weights,strain_tolerance,strain_var_mat_previousinc,IsProj,RoutineID);
        [~, Res_F_E, Res_F_F] = func_partitionf(fixnodes_applied,Res_F);
        Res_F_E_rct=Res_F_E+f_rct_disp_ebc;
        
        if SolverID == 2  % For non-local gradient only
            % Partition the residual vector (free) based on the contribution of
            % the displacement and non-local strains 
            [Res_F_F_u_f,Res_F_F_e_star] = func_unroll_free_dofs_UAL_gradient(ID_dofs_list_u_p,ID_dofs_list_u_f,ID_dofs_list_nl_strain,Res_F_F,Res_F_E);     
            % Paritition Consistent Tangent matrix J
            [Jpf,Jpp,Jff,Jfp,Jue_p,Jue_f,Jeu_p,Jeu_f,Jee] = func_sort_stiffness_UAL_gradient(J,ID_dofs_list_u_p,ID_dofs_list_u_f,ID_dofs_list_nl_strain,ID_dofs_list_disp,fixnodes_applied);
            % Calculate Matrix 1 (if applicable)
            [~,struct_sub_Matrix_1] = func_calc_Matrix_1_UAL_gradient(Jff,Jpp,Jpf,Jfp,Jue_p,Jue_f,Jeu_p,Jeu_f,Jee,Applied_Displacement_Load,Beta,delta_u_bar_current,delta_u_f_current,delta_f_rct_disp_ebc_current,M_bar_u,M_u,M_bar_e,M_e);
        
        end

        % Calculate and store residual norm
        [g_constraint] = func_ArcLength_constraint_UAL(SolverID, delta_e_nl_current, delta_u_bar_current,delta_u_f_current,delta_f_rct_disp_ebc_current,ArcLength,Beta);        
        residual=[Res_F_E_rct;Res_F_F;g_constraint];
        residual_norm(iteration)=norm(residual,2);      

        % If there isn't a significant reduction in Res_norm within 10
        % iterations, go back and reduce the ArcLength
        if ((iteration>11) && ((residual_norm(iteration-10)/residual_norm(iteration))<10)) || residual_norm(iteration)>1e6
            Restore_arc_length    = ArcLength;  
            Restore_residual_norm = residual_norm;
            cd (save_path)
            format shortg
            last_saved_loadstep=sprintf('Increment = %d, Load percentage = %.3f .mat',increment-1,Load_percentage_last_saved);
            load (last_saved_loadstep)
            cd (main_file_path)
            ArcLength     = Restore_arc_length;   
            residual_norm = Restore_residual_norm;
            return
        end

        % If the residual error is less than the convergence threshold then the
        % solution is converged, otherwise proceed to the next iteration
        if (residual_norm(iteration) < tolerance) || (iteration == max_accept_iter)                               
            if (residual_norm(iteration) < tolerance)
                % Calcualte Reactions based on the geometry of your domain and boundary conditions
                f_ext_E=f_rct_disp_ebc;
                if  model_name=='SSNT_Coarse' || model_name=='SSNT_Fine' || model_name=='SNS' || model_name=='TNT' 
                    Reactions_x = sum(abs(f_ext_E(fixnodes_applied(2,:) == 1 & fixnodes_applied(3,:) == 0)));
                    Reactions_y = sum(abs(f_ext_E(fixnodes_applied(2,:) == 2 & fixnodes_applied(3,:) == 0)));
                elseif model_name=='SNT'
                    Reactions_x = sum(abs(f_ext_E(fixnodes_applied(2,:) == 1 & fixnodes_applied(3,:) ~= 0)));
                    Reactions_y = sum(abs(f_ext_E(fixnodes_applied(2,:) == 2 & fixnodes_applied(3,:) ~= 0)));
                end

                % Check the force equilibrium
                Sum_Fx = sum(f_ext_E(fixnodes_applied(2,:) == 1))
                Sum_Fy = sum(f_ext_E(fixnodes_applied(2,:) == 2))

                disp("End of load increment no.: " + num2str(increment))
                disp("iteration " + num2str(iteration) + " Residual norm: " + "  " + residual_norm(iteration))
                disp("========================================")

                % Save converged deltas
                delta_m_bar_conv            = delta_m_bar;
                delta_u_bar_conv            = delta_u_bar;
                delta_u_f_conv              = delta_u_f;
                delta_f_rct_disp_ebc_conv   = delta_f_rct_disp_ebc;
                
                if SolverID == 2 
                    delta_e_nl_conv         = delta_e_nl;
                end

                % Save converged state variables
                m_bar_conv             = m_bar;
                u_bar_conv             = u_bar;
                u_f_conv               = u_f;
                f_rct_disp_ebc_conv    = f_rct_disp_ebc;
                if SolverID == 2 
                    e_nl_conv          = e_nl;
                end

                % Save converged dofs 
                dofs_stored=dofs;

                % Saved history variables 
                history_var_mat_conv = history_var_mat;
                strain_var_mat_conv  = strain_var_mat;
                convergance_flag     = 1;
            else
                Sum_Fx = [];
                Sum_Fy = [];
                Reactions_x = [];
                Reactions_y = [];
                convergance_flag     = 0;
            end

            break
        end

    end

elseif IsProj == 1
    % ---------------------------------------------------------------------
    %%%%%%%%%%%%%%%%%%%%% PLOTTING PROPERTIES ROUTINE %%%%%%%%%%%%%%%%%%%%%
    % ---------------------------------------------------------------------

    % Create empty entries for solver variables
    [convergance_flag,delta_e_nl_conv,e_nl_conv,delta_m_bar,increment,J, dofs_stored, residual_norm, g_constraint,Res_F_E_rct,Res_F_F,iteration, history_var_mat_conv, strain_var_mat_conv, gausspoints_prop_mat, nodes_prop_mat,Reactions_x, Reactions_y,delta_m_bar_conv,delta_u_bar_conv,delta_u_f_conv,delta_f_rct_disp_ebc_conv,m_bar_conv,u_bar_conv,u_f_conv,f_rct_disp_ebc_conv,ArcLength,Load_percentage_last_saved] = deal([]);
       
    % Compute element/nodal properties
    [~, ~, ~, ~, ~, ~, ~, gausspoints_prop_mat, nodes_prop_mat,~] = func_globalstiffness(model_name,dofs,Delastic,history_var_mat_previousinc,num_elem_at_node,n_hood,weights,strain_tolerance,strain_var_mat_previousinc,IsProj,RoutineID);

else

    disp("Check your IsProj variable")

end




end


