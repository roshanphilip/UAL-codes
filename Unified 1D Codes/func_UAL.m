function [delta_m_bar_con,delta_u_p_con,delta_u_f_con,delta_e_nl_con, delta_f_rct_con,m_bar_con, u_p_con, u_f_con, e_nl_con, f_rct_con, n, k, storage_damage_loadstep, storage_strain_loadstep, storage_nonlocal_strain_loadstep, residual_norm, ArcLength] = func_UAL(ST, delta_m_bar_0, constraint_type, n_dofs, c, delta_m_bar_con,delta_u_p_con,delta_u_f_con,delta_e_nl_con, delta_f_rct_con, m_bar_con, u_p_con, u_f_con, e_nl_con, f_rct_con,n_NodesPerElement,weight_integrationpoint,a,b,n_TotalElements,Applied_Displacement_Load,tolerance,k_max,n,ArcLength,DamageThreshholdStrain,n_TotalNodes,n_IntegrationPoints,ConstitutiveMatrix,YoungsModulus,Globalcoordinates_EachNode,connect_LocaltoGlobal_Nodenumbers,n_Global,storage_damage_loadstep_con,ID_prescribed_u,ID_free_u,ID_u,ID_nl_strain,ID_free,SolverID)
% This function contains the UAL routine

M_bar    = size(ID_prescribed_u,1);
M        = size(ID_free_u,1);

% Recover last converged state variable values
m_bar    =   m_bar_con;
u_p      =   u_p_con;
u_f      =   u_f_con;
f_rct    =   f_rct_con;
e_nl     =   e_nl_con;

% Initialize iteration counter
k        =   0;

% Set residual norm to a high value to start the iterations
residual_norm=2*tolerance;

while residual_norm(1,end)>tolerance && k<k_max
    % Update iteration counter
    k=k+1;

    % Compute Beta
    if SolverID == 1     % Local damage
        dofs(ID_prescribed_u,1) = u_p;
        dofs(ID_free_u,1)       = u_f;
    elseif SolverID ==2  % Gradient damage
        dofs(ID_prescribed_u,1) = u_p;
        dofs(ID_free_u,1)       = u_f;
        dofs(ID_nl_strain,1)    = e_nl;
    end
    [Beta] = func_calculate_Beta(ST, constraint_type,M_bar,n_dofs, c, n_TotalElements,n_TotalNodes,n_NodesPerElement,connect_LocaltoGlobal_Nodenumbers,Globalcoordinates_EachNode,n_IntegrationPoints,ConstitutiveMatrix,weight_integrationpoint,DamageThreshholdStrain,a,b,n_Global,storage_damage_loadstep_con,YoungsModulus,dofs,ID_u,ID_nl_strain,ID_prescribed_u,SolverID);

    % A. Only for the 1st iteration
    if k==1

        % 1. Only for the 1st loadstep
        if n==1

            % Calculate consistent tangent matrix
            if SolverID == 1     % Local damage
                dofs(ID_prescribed_u,1) = u_p;
                dofs(ID_free_u,1)       = u_f;
            elseif SolverID ==2  % Gradient damage
                dofs(ID_prescribed_u,1) = u_p;
                dofs(ID_free_u,1)       = u_f;
                dofs(ID_nl_strain,1)    = e_nl;
            end
            [~,J,~,~,~,~,~,~] = func_CalculateGlobalMatrices(ST, n_dofs, c, n_TotalElements,n_TotalNodes,n_NodesPerElement,connect_LocaltoGlobal_Nodenumbers,Globalcoordinates_EachNode,n_IntegrationPoints,ConstitutiveMatrix,weight_integrationpoint,DamageThreshholdStrain,a,b,n_Global,storage_damage_loadstep_con,YoungsModulus,dofs,ID_u,ID_nl_strain,SolverID);

            % Calculate predictor values
            [delta_m_bar,delta_u_p,delta_u_f,delta_f_rct,delta_e_nl] = func_delta_predictors_UAL(n, Applied_Displacement_Load,delta_m_bar_0,delta_m_bar_con,J,SolverID,ID_prescribed_u,ID_free_u,ID_nl_strain);

        else
            % Recover last converged deltas
            delta_m_bar                   =   delta_m_bar_con;
            delta_u_p                     =   delta_u_p_con;
            delta_u_f                     =   delta_u_f_con;
            delta_f_rct                   =   delta_f_rct_con;
            delta_e_nl                    =   delta_e_nl_con;
        end

    else

        if SolverID == 1
            % Calculate corrector values
            [del_f_rct, del_u_f, del_m_bar, del_e_nl] = func_calc_dels_PC_scheme(Residual_Global, g, delta_u_p, delta_u_f, delta_f_rct, J, Beta, Applied_Displacement_Load, ID_prescribed_u, ID_free_u, f_rct);
        elseif SolverID == 2
            [del_f_rct, del_u_f, del_m_bar, del_e_nl] = func_solve_system_partitioned_gradient(struct_sub_Matrix_1,Residual_Global, g, M_bar, M, ID_prescribed_u, ID_free_u, ID_u, ID_nl_strain, ID_free, f_rct);
        end
        % Update deltas
        delta_m_bar                       = delta_m_bar    + del_m_bar;
        delta_u_p                         = delta_m_bar    * Applied_Displacement_Load;
        delta_u_f                         = delta_u_f      + del_u_f;
        delta_f_rct                       = delta_f_rct    + del_f_rct;
        delta_e_nl                        = delta_e_nl     + del_e_nl;
    end

    % Save current state variables
    delta_m_bar_current                   =   delta_m_bar;
    delta_u_p_current                     =   delta_u_p;
    delta_u_f_current                     =   delta_u_f;
    delta_f_rct_current                   =   delta_f_rct;
    delta_e_nl_current                    =   delta_e_nl;

    % Update state variables
    m_bar                                 =   m_bar_con + delta_m_bar_current;
    u_p                                   =   m_bar     * Applied_Displacement_Load;
    u_f                                   =   u_f_con   + delta_u_f_current;
    f_rct                                 =   f_rct_con + delta_f_rct_current;
    e_nl                                  =   e_nl_con  + delta_e_nl_current;

    % Calculate residuals
    if SolverID == 1     % Local damage
        dofs(ID_prescribed_u,1) = u_p;
        dofs(ID_free_u,1)       = u_f;
    elseif SolverID ==2  % Gradient damage
        dofs(ID_prescribed_u,1) = u_p;
        dofs(ID_free_u,1)       = u_f;
        dofs(ID_nl_strain,1)    = e_nl;
    end
    [~,J,~,~,Residual_Global,storage_damage_loadstep,storage_strain_loadstep,storage_nonlocal_strain_loadstep] = func_CalculateGlobalMatrices(ST, n_dofs, c, n_TotalElements,n_TotalNodes,n_NodesPerElement,connect_LocaltoGlobal_Nodenumbers,Globalcoordinates_EachNode,n_IntegrationPoints,ConstitutiveMatrix,weight_integrationpoint,DamageThreshholdStrain,a,b,n_Global,storage_damage_loadstep_con,YoungsModulus,dofs,ID_u,ID_nl_strain,SolverID);
    [g] = func_UAL_arclength_eq(delta_u_p_current,delta_u_f_current,delta_f_rct_current,delta_e_nl_current,ArcLength,Beta,SolverID);
    Res_u_p  = Residual_Global(ID_prescribed_u,1) + f_rct;
    Res_f  = Residual_Global(ID_free,1);
    residual = [Res_u_p;Res_f;g];

    if SolverID == 2
        [struct_sub_Matrix_1] = func_calc_Matrix_1(J,Applied_Displacement_Load,Beta,delta_u_p_current,delta_u_f_current,delta_f_rct_current,M_bar,M,ID_prescribed_u,ID_free_u,ID_nl_strain);
    end

    residual_norm(1,k) = norm(residual,2);

    if residual_norm(1,k) > 1e2
        break
    end

end

% Check convergence
if residual_norm(1,end)<tolerance
    % Print converged residual value
    residual_norm(1,end);

    % Save state variables and deltas
    delta_m_bar_con     =   delta_m_bar_current;
    delta_u_p_con       =   delta_u_p_current;
    delta_u_f_con       =   delta_u_f_current;
    delta_f_rct_con     =   delta_f_rct_current;
    delta_e_nl_con      =   delta_e_nl_current;

    m_bar_con           =   m_bar;
    u_p_con             =   u_p;
    u_f_con             =   u_f;
    f_rct_con           =   f_rct;
    e_nl_con            =   e_nl;
end


end