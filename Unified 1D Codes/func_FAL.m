function [delta_dofs_con, delta_lambda_con,delta_u_p_con,delta_u_f_con,delta_e_nl_con, delta_f_rct_con, lambda_con, u_p_con, u_f_con, e_nl_con, f_reaction, n, k, storage_damage_loadstep, storage_strain_loadstep, storage_nonlocal_strain_loadstep, residual_norm, ArcLength] = func_FAL(ST, q, n_FreeNodes, n_EssentialNodes, constraint_type, n_dofs, c, delta_lambda_con,delta_u_p_con,delta_u_f_con,delta_e_nl_con, delta_f_rct_con, lambda_con, u_p_con, u_f_con, e_nl_con, f_reaction,n_NodesPerElement,weight_integrationpoint,a,b,n_TotalElements,tolerance,k_max,n,ArcLength,DamageThreshholdStrain,n_TotalNodes,n_IntegrationPoints,ConstitutiveMatrix,YoungsModulus,Globalcoordinates_EachNode,connect_LocaltoGlobal_Nodenumbers,n_Global,storage_damage_loadstep_con,ID_prescribed_u,ID_free_u,ID_u,ID_nl_strain,ID_free,SolverID,delta_dofs_con)
% This function contains the FAL routine

% Set k = 1
k = 1;

% Set delta_(.) = 0
delta_dofs                    = zeros(size(ID_prescribed_u,1)+size(ID_free,1),1);
delta_lambda                  = 0;

% Set del_u_A = 0
del_dofs_A_free               = zeros(size(ID_free,1),1);

% Calculate del_u_B
q_e   = q(ID_prescribed_u,1);
q_f   = q(ID_free,1);

if SolverID == 1       % Local damage
    dofs(ID_prescribed_u,1)  = u_p_con;
    dofs(ID_free_u,1)  = u_f_con;
elseif SolverID == 2   % Gradient damage
    dofs(ID_prescribed_u,1)  = u_p_con;
    dofs(ID_free_u,1)  = u_f_con;
    dofs(ID_nl_strain,1) = e_nl_con;
end


[~, J, ~, ~, ~, ~, ~, ~] = func_CalculateGlobalMatrices(ST, n_dofs, c, n_TotalElements,n_TotalNodes,n_NodesPerElement,connect_LocaltoGlobal_Nodenumbers,Globalcoordinates_EachNode,n_IntegrationPoints,ConstitutiveMatrix,weight_integrationpoint,DamageThreshholdStrain,a,b,n_Global,storage_damage_loadstep_con,YoungsModulus,dofs,ID_u,ID_nl_strain,SolverID);
Jff = J(ID_free,ID_free);

del_dofs_B_free=Jff\q_f;

% Calculate Beta
[Beta] = func_calculate_Beta(ST, constraint_type,size(ID_prescribed_u,1),n_dofs, c, n_TotalElements,n_TotalNodes,n_NodesPerElement,connect_LocaltoGlobal_Nodenumbers,Globalcoordinates_EachNode,n_IntegrationPoints,ConstitutiveMatrix,weight_integrationpoint,DamageThreshholdStrain,a,b,n_Global,storage_damage_loadstep_con,YoungsModulus,dofs,ID_u,ID_nl_strain,ID_prescribed_u,SolverID);

% Calculate del_lambda
[del_lambda,imaginary_roots] = func_correctors_FAL(Jff,delta_lambda_con,delta_u_p_con,delta_dofs_con(ID_free,1),delta_e_nl_con,1,del_dofs_A_free,del_dofs_B_free,q_f,Beta,delta_lambda,delta_dofs,ArcLength,1,ID_prescribed_u,ID_free);

% If imaginary roots are present, go back and decrease the ArcLength
if imaginary_roots == 1
    residual_norm = 2 * tolerance;
    return
end

% Calculate del_dofs
del_dofs_free               = del_dofs_A_free + (del_dofs_B_free * del_lambda);
del_dofs(ID_prescribed_u,1) = 0;
del_dofs(ID_free,1)         = del_dofs_free;

% Update state variables
dofs_con(ID_prescribed_u,1) = u_p_con;
dofs_con(ID_free_u,1)       = u_f_con;
u                           = dofs_con(ID_u,1)     + del_dofs(ID_u,1);
dofs(ID_u,1)                = u;
if SolverID == 2 
    e_nl                 = e_nl_con  + del_dofs(ID_nl_strain,1);
    dofs(ID_nl_strain,1) = e_nl;
end
lambda                   = lambda_con + del_lambda;

% Calculate Residual
[K, J, F_internal_Global, ~, ~, storage_damage_loadstep, storage_strain_loadstep, storage_nonlocal_strain_loadstep] = func_CalculateGlobalMatrices(ST, n_dofs, c, n_TotalElements,n_TotalNodes,n_NodesPerElement,connect_LocaltoGlobal_Nodenumbers,Globalcoordinates_EachNode,n_IntegrationPoints,ConstitutiveMatrix,weight_integrationpoint,DamageThreshholdStrain,a,b,n_Global,storage_damage_loadstep_con,YoungsModulus,dofs,ID_u,ID_nl_strain,SolverID);
F_int_f                                  = F_internal_Global(ID_free,1);
[gamma]                                  = func_FAL_arclength_eq(Beta,del_dofs,del_lambda,delta_dofs,delta_lambda,q,ArcLength);
Residual_free                            = F_int_f - (lambda * q_f);
Residual_u(ID_prescribed_u,1)            = 0;
Residual_u(ID_free,1)                    = Residual_free;
Residual                                 = [Residual_u;gamma];
residual_norm(1,k)                       = norm(Residual);

% Check convergence
if residual_norm(1,end) <= tolerance
    % Save converged variables
    u_p_con          = dofs(ID_prescribed_u,1);
    u_f_con          = dofs(ID_free_u,1);
    dofs_con         = dofs;
    lambda_con       = lambda;
    delta_u_p_con    = delta_dofs(ID_prescribed_u,1);
    delta_u_f_con    = delta_dofs(ID_free_u,1);
    delta_lambda_con = delta_lambda;
    delta_dofs_con   = delta_dofs;
    if SolverID == 2
        e_nl_con       = dofs(ID_nl_strain,1);
        delta_e_nl_con = delta_dofs(ID_nl_strain,1);
    end

    % Calculate reaction force
    f_reaction       = F_internal_Global(ID_u(end,1),1);
    f_reaction2      = q * lambda;

    % Print converged residual value
    residual_norm(1,end);

end

while residual_norm(1,end) > tolerance && k<k_max

    % Update iteration counter
    k = k + 1;

    % Update deltas
    delta_dofs        = del_dofs;
    delta_lambda      = del_lambda;

    [~, J, F_internal_Global, ~, ~, storage_damage_loadstep, storage_strain_loadstep, storage_nonlocal_strain_loadstep] = func_CalculateGlobalMatrices(ST, n_dofs, c, n_TotalElements,n_TotalNodes,n_NodesPerElement,connect_LocaltoGlobal_Nodenumbers,Globalcoordinates_EachNode,n_IntegrationPoints,ConstitutiveMatrix,weight_integrationpoint,DamageThreshholdStrain,a,b,n_Global,storage_damage_loadstep_con,YoungsModulus,dofs,ID_u,ID_nl_strain,SolverID);

    % Calculate del_dofs_A_free
    Jff            = J(ID_free,ID_free);
    F_int_f        = F_internal_Global(ID_free,1);
    del_dofs_A_free   = -Jff\(F_int_f-((lambda)*q_f));

    % Calculate del_dofs_B_free
    del_dofs_B_free   = Jff\q_f;

    % Calculate del_lambda
    [del_lambda,imaginary_roots] = func_correctors_FAL(Jff,delta_lambda_con,delta_u_p_con,delta_dofs_con(ID_free,1),delta_e_nl_con,1,del_dofs_A_free,del_dofs_B_free,q_f,Beta,delta_lambda,delta_dofs,ArcLength,1,ID_prescribed_u,ID_free);

    % If imaginary roots are present, go back and decrease the ArcLength
    if imaginary_roots == 1
        residual_norm(1,k) = 2 * tolerance;
        return
    else

        % Calculate del_dofs
        del_dofs_free               = del_dofs_A_free + (del_dofs_B_free * del_lambda);
        del_dofs(ID_prescribed_u,1) = 0;
        del_dofs(ID_free,1)         = del_dofs_free;

        % Update state variables
        dofs_con(ID_prescribed_u,1) = u_p_con           +  delta_dofs(ID_prescribed_u,1)  + del_dofs(ID_prescribed_u,1);
        dofs_con(ID_free_u,1)       = u_f_con           +  delta_dofs(ID_free_u,1)        + del_dofs(ID_free_u,1);
        u                           = dofs_con(ID_u,1);
        dofs(ID_u,1)                = u;
        if SolverID == 2
            e_nl                 = e_nl_con  + delta_dofs(ID_nl_strain,1) + del_dofs(ID_nl_strain,1);
            dofs(ID_nl_strain,1) = e_nl;
        end
        lambda                   = lambda_con + delta_lambda + del_lambda;

        % Calculate Residual
        [K, J, F_internal_Global, ~, ~, storage_damage_loadstep, storage_strain_loadstep, storage_nonlocal_strain_loadstep] = func_CalculateGlobalMatrices(ST, n_dofs, c, n_TotalElements,n_TotalNodes,n_NodesPerElement,connect_LocaltoGlobal_Nodenumbers,Globalcoordinates_EachNode,n_IntegrationPoints,ConstitutiveMatrix,weight_integrationpoint,DamageThreshholdStrain,a,b,n_Global,storage_damage_loadstep_con,YoungsModulus,dofs,ID_u,ID_nl_strain,SolverID);
        F_int_f                                  = F_internal_Global(ID_free,1);
        [gamma]                                  = func_FAL_arclength_eq(Beta,del_dofs,del_lambda,delta_dofs,delta_lambda,q,ArcLength);
        Residual_free                            = F_int_f - (lambda * q_f);
        Residual_u(ID_prescribed_u,1)            = 0;
        Residual_u(ID_free,1)                    = Residual_free;
        Residual                                 = [Residual_u;gamma];
        residual_norm(1,k)                       = norm(Residual);


        % Check convergence
        if residual_norm(1,end) <= tolerance
            % Save converged variables
            u_p_con          = dofs(ID_prescribed_u,1);
            u_f_con          = dofs(ID_free_u,1);
            dofs_con         = dofs;
            lambda_con       = lambda;
            delta_u_p_con    = delta_dofs(ID_prescribed_u,1);
            delta_u_f_con    = delta_dofs(ID_free_u,1);
            delta_lambda_con = delta_lambda;
            delta_dofs_con   = delta_dofs;
            if SolverID == 2
                e_nl_con       = dofs(ID_nl_strain,1);
                delta_e_nl_con = delta_dofs(ID_nl_strain,1);
            end
            
            % Calculate reaction force
            f_reaction       = F_internal_Global(ID_u(end,1),1);
            f_reaction2      = q * lambda;

            % Print converged residual value
            residual_norm(1,end);

        end
    end
end
end