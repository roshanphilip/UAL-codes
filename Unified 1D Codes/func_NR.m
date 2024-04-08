function [dofs_con,f_reaction,n,k,storage_damage_loadstep,storage_strain_loadstep,storage_nonlocal_strain_loadstep,residual_norm,delta_lf] = func_NR(ST, n_dofs, c, dofs_con, lf,n_NodesPerElement,weight_integrationpoint,a,b,n_TotalElements,Applied_Displacement_Load,tolerance,k_max,n,delta_lf,DamageThreshholdStrain,n_TotalNodes,n_IntegrationPoints,f_reaction,ConstitutiveMatrix,YoungsModulus,Globalcoordinates_EachNode,connect_LocaltoGlobal_Nodenumbers,n_Global,storage_damage_loadstep_con,ID_prescribed_u,ID_free_u,ID_u,ID_nl_strain,ID_free,SolverID)
% This function contains the NR routine

% Reset iteration counter
k          = 1;

% Update loadfactor
lf         = lf + delta_lf;

% Update essential nodes
del_dofs_p              = lf * Applied_Displacement_Load;
dofs                    = dofs_con;
dofs(ID_prescribed_u,1) = dofs(ID_prescribed_u,1) + del_dofs_p;

% Calculate consistent tangent matrix 
[K,J,~,~,Residual_Global,storage_damage_loadstep,storage_strain_loadstep,storage_nonlocal_strain_loadstep] = func_CalculateGlobalMatrices(ST, n_dofs, c, n_TotalElements,n_TotalNodes,n_NodesPerElement,connect_LocaltoGlobal_Nodenumbers,Globalcoordinates_EachNode,n_IntegrationPoints,ConstitutiveMatrix,weight_integrationpoint,DamageThreshholdStrain,a,b,n_Global,storage_damage_loadstep_con,YoungsModulus,dofs,ID_u,ID_nl_strain,SolverID);

% Sort Global matrices
Jff   = J(ID_free,ID_free);
Res_f = Residual_Global(ID_free,1);

% Calculate residuals
residual_norm(1,k) = norm(Res_f,2);

while residual_norm(1,end)>tolerance && k<k_max

    % Update iteration counter
    k = k+1;

    % Update position of free nodes
    del_dofs_f      = -Jff\Res_f;    
    dofs(ID_free,1) = dofs(ID_free,1) + del_dofs_f;    

    % Calculate consistent tangent matrix 
    [K,J,~,~,Residual_Global,storage_damage_loadstep,storage_strain_loadstep,storage_nonlocal_strain_loadstep] = func_CalculateGlobalMatrices(ST, n_dofs, c, n_TotalElements,n_TotalNodes,n_NodesPerElement,connect_LocaltoGlobal_Nodenumbers,Globalcoordinates_EachNode,n_IntegrationPoints,ConstitutiveMatrix,weight_integrationpoint,DamageThreshholdStrain,a,b,n_Global,storage_damage_loadstep_con,YoungsModulus,dofs,ID_u,ID_nl_strain,SolverID);

    % Sort Global matrices
    Jff = J(ID_free,ID_free);
    Res_f = Residual_Global(ID_free,1);

    % Calculate residuals
    residual_norm(1,k) = norm(Res_f,2);

    % Exit increment if residuals are too large
    if residual_norm(1,k)>1e3
        break
    end

end

% Check convergence
if residual_norm(1,end)<tolerance
    % Print converged residual value
    residual_norm(1,end);

    % Paritition stiffness matrix
    Kee = K(ID_prescribed_u,ID_prescribed_u);
    Kef = K(ID_prescribed_u,ID_free_u);

    % Save displacement and reactions
    dofs_con           = dofs;
    f_reaction         = Kee * dofs_con(ID_prescribed_u,1) + Kef * dofs_con(ID_free_u,1);

    clc
else
    f_reaction         = [];
end

end