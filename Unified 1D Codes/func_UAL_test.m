function [f_reaction, f_reaction_essential_con,m_bar_con,dofs_con,n,k,storage_damage_loadstep,storage_strain_loadstep,storage_nonlocal_strain_loadstep,residual_norm,delta_u_bar_con,delta_u_f_con,delta_f_reaction_essential_con,delta_m_bar_con,delta_e_nl_con,ArcLength] = func_UAL_test(f_reaction_essential_con,m_bar_con,ArcLength,delta_u_bar_con,delta_u_f_con,delta_f_reaction_essential_con,delta_m_bar_con,delta_e_nl_con, constraint_type, n_dofs, c, dofs_con, n_NodesPerElement, weight_integrationpoint, a, b, n_TotalElements, Applied_Displacement_Load, tolerance, k_max, n, DamageThreshholdStrain, n_TotalNodes, n_IntegrationPoints, ConstitutiveMatrix, YoungsModulus, Globalcoordinates_EachNode, connect_LocaltoGlobal_Nodenumbers, n_Global, storage_damage_loadstep_con, ID_prescribed_u, ID_free_u, ID_u,ID_nl_strain, ID_free, SolverID, delta_m_bar_0)
% This function uses the NR solver to find the solution

if SolverID == 1
    [e_nl, e_nl_con, delta_e_nl, delta_e_nl_current, delta_e_nl_con] = deal(zeros());
end

% Recover last converged values of dofs
dofs                  = dofs_con;
m_bar                 = m_bar_con;
f_reaction_essential  = f_reaction_essential_con;

% Unroll the u_f_con vector
u_f_con               = dofs_con(ID_free_u,1);

if SolverID == 2
    % Unroll the e_nl_con vector
    e_nl_con              = dofs_con(ID_nl_strain,1);
end

% Update dofs vector
u_bar                 = dofs(ID_prescribed_u);
u_f                   = dofs(ID_free_u);
if SolverID == 2
    e_nl                  = dofs(ID_nl_strain);
end

%-------------------------------------------------------------------------%

% Reset iteration counter
k             = 0;

% Set residual norm to a high value to start the iterations
residual_norm = 2 * tolerance;

% Iteration loop
while residual_norm(1,end)>tolerance && k<k_max

    % Update iteration counter
    k=k+1;

    % Compute Beta
    [Beta] = func_calculate_Beta(constraint_type,size(ID_prescribed_u,1),n_dofs, c, n_TotalElements,n_TotalNodes,n_NodesPerElement,connect_LocaltoGlobal_Nodenumbers,Globalcoordinates_EachNode,n_IntegrationPoints,ConstitutiveMatrix,weight_integrationpoint,DamageThreshholdStrain,a,b,n_Global,storage_damage_loadstep_con,YoungsModulus,dofs,ID_u,ID_nl_strain,ID_prescribed_u,SolverID);

    % A. Only for the 1st iteration
    if k==1
        % 1. Only for the 1st loadstep
        if n==1
            % Calculate consistent tangent matrix
            [~,J,~,~,~,~,~,~] = func_CalculateGlobalMatrices(n_dofs, c, n_TotalElements,n_TotalNodes,n_NodesPerElement,connect_LocaltoGlobal_Nodenumbers,Globalcoordinates_EachNode,n_IntegrationPoints,ConstitutiveMatrix,weight_integrationpoint,DamageThreshholdStrain,a,b,n_Global,storage_damage_loadstep_con,YoungsModulus,dofs,ID_u,ID_nl_strain,SolverID);

            % Calculate predictor values for deltas
            [delta_m_bar,delta_u_bar,delta_u_f,delta_f_reaction_essential,delta_e_nl] = func_delta_predictors_UAL(n, Applied_Displacement_Load,delta_m_bar_0, delta_m_bar_con, J,SolverID,ID_prescribed_u,ID_free_u,ID_nl_strain);       
        else

            % Calculate consistent tangent matrix
            [~,J,~,~,~,~,~,~] = func_CalculateGlobalMatrices(n_dofs, c, n_TotalElements,n_TotalNodes,n_NodesPerElement,connect_LocaltoGlobal_Nodenumbers,Globalcoordinates_EachNode,n_IntegrationPoints,ConstitutiveMatrix,weight_integrationpoint,DamageThreshholdStrain,a,b,n_Global,storage_damage_loadstep_con,YoungsModulus,dofs,ID_u,ID_nl_strain,SolverID);

            % Calculate predictor values for deltas
            [delta_m_bar,delta_u_bar,delta_u_f,delta_f_reaction_essential,delta_e_nl] = func_delta_predictors_UAL(n, Applied_Displacement_Load,delta_m_bar_0, delta_m_bar_con, J,SolverID,ID_prescribed_u,ID_free_u,ID_nl_strain);       
       

            % % Recover last converged deltas
            % delta_m_bar                 =  delta_m_bar_con;
            % delta_u_bar                 =  delta_u_bar_con;
            % delta_u_f                   =  delta_u_f_con;
            % delta_f_reaction_essential  =  delta_f_reaction_essential_con;
            % if SolverID == 2
            %     delta_e_nl              =  delta_e_nl_con;
            % end
        end

        % B. Any iteration after the 1st
    else
        % Calculate dels
        if SolverID == 1        % Local damage
            del_e_nl                                      = [];
            [del_f_reaction, del_u_f, del_m_bar]          = func_solve_system_equations_consistent_partitioned(Residual_Global, g, delta_u_bar, delta_u_f, delta_f_reaction_essential, J, Beta, Applied_Displacement_Load, ID_prescribed_u, ID_free_u);
            
            if isnan(del_m_bar)  == 1
                break
            end

        elseif SolverID == 2    % Gradient damage
            [del_f_reaction, del_u_f, del_m_bar,del_e_nl] = func_solve_system_partitioned_gradient(struct_sub_Matrix_1, Residual_Global, g, size(ID_prescribed_u,1), size(ID_free_u,1), ID_prescribed_u, ID_free_u, ID_u, ID_nl_strain, ID_free);
        
            if isnan(del_m_bar)  == 1
                break
            end

        end

        % Update deltas
        delta_m_bar                 =  delta_m_bar                 + del_m_bar;                        % size: [1,1]
        delta_u_bar                 =  delta_m_bar                 * Applied_Displacement_Load;        % size: [M_bar,1]
        delta_u_f                   =  delta_u_f                   + del_u_f;                          % size: [M,1]
        delta_f_reaction_essential  =  delta_f_reaction_essential  + del_f_reaction;                   % size: [M_bar,1]
        if SolverID == 2
            delta_e_nl              =  delta_e_nl                  + del_e_nl;                         % size: [M_bar+M,1]
        end
    end

    % Save the current delta values
    delta_m_bar_current                 =  delta_m_bar;                                                % size: [1,1] 
    delta_u_bar_current                 =  delta_m_bar_current * Applied_Displacement_Load;            % size: [M_bar,1]
    delta_u_f_current                   =  delta_u_f;                                                  % size: [M,1]
    delta_f_reaction_essential_current  =  delta_f_reaction_essential;                                 % size: [M_bar,1]
    if SolverID == 2
        delta_e_nl_current              =  delta_e_nl;                                                 % size: [M_bar+M,1]                
    end

    % Update state variables
    m_bar                 =  m_bar_con                + delta_m_bar_current;                           % size: [1,1]
    u_bar                 =  m_bar                    * Applied_Displacement_Load;                     % size: [M_bar,1]
    u_f                   =  u_f_con                  + delta_u_f_current;                             % size: [M,1]
    f_reaction_essential  =  f_reaction_essential_con + delta_f_reaction_essential_current;            % size: [M_bar,1]
    if SolverID == 2
        e_nl              =  e_nl_con                 + delta_e_nl_current;                            % size: [M_bar+M,1]
    end

    % Update dofs vector
    dofs(ID_prescribed_u)  =  u_bar;
    dofs(ID_free_u)        =  u_f;
    if SolverID == 2
        dofs(ID_nl_strain) =  e_nl;
    end

    % Calculate residuals
    % Calculate g
    [g] = func_UAL_arclength_eq(delta_u_bar_current,delta_u_f_current,delta_f_reaction_essential_current,delta_e_nl_current,ArcLength,Beta,SolverID);
    % Calculate consistent tangent matrix
    [~,J,~,~,Residual_Global,storage_damage_loadstep,storage_strain_loadstep,storage_nonlocal_strain_loadstep] = func_CalculateGlobalMatrices(n_dofs, c, n_TotalElements,n_TotalNodes,n_NodesPerElement,connect_LocaltoGlobal_Nodenumbers,Globalcoordinates_EachNode,n_IntegrationPoints,ConstitutiveMatrix,weight_integrationpoint,DamageThreshholdStrain,a,b,n_Global,storage_damage_loadstep_con,YoungsModulus,dofs,ID_u,ID_nl_strain,SolverID);

    if SolverID == 2 % Gradient damage
        % Calculate Matrix 1 (if applicable)
        [struct_sub_Matrix_1]       = func_calc_Matrix_1(J, Applied_Displacement_Load,Beta,delta_u_bar_current,delta_u_f_current,delta_f_reaction_essential_current,size(ID_prescribed_u,1),size(ID_free,1),ID_prescribed_u,ID_free_u,ID_nl_strain);
    end

    % Calculate residual norm
    residual_norm(1,k)              = norm([(Residual_Global(ID_prescribed_u,1)+f_reaction_essential);Residual_Global(ID_free,1);g]);    

end

% Check convergence
if residual_norm(1,end)<tolerance && k<k_max
    % Print converged residual value
    residual_norm(1,end)

    % Save converged displacement and reactions
    dofs_con                       =   dofs;
    m_bar_con                      =   m_bar;
    f_reaction_essential_con       =   f_reaction_essential;

    % Save converged deltas
    delta_m_bar_con                =   delta_m_bar;
    delta_u_bar_con                =   delta_u_bar;
    delta_u_f_con                  =   delta_u_f;
    delta_f_reaction_essential_con =   delta_f_reaction_essential;

    f_reaction                     =   f_reaction_essential;
else
    f_reaction                     =   []; 
end
end