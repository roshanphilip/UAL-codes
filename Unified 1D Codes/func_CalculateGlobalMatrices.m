function [K, J, F_internal_Global, F_external_Global, Residual_Global, storage_damage_loadstep, storage_strain_loadstep, storage_nonlocal_strain_loadstep] = func_CalculateGlobalMatrices(ST, n_dofs, c, n_TotalElements,n_TotalNodes,n_NodesPerElement,connect_LocaltoGlobal_Nodenumbers,Globalcoordinates_EachNode,n_IntegrationPoints,ConstitutiveMatrix,weight_integrationpoint,DamageThreshholdStrain,a,b,n_Global,storage_damage_loadstep_con,YoungsModulus,dofs,ID_u,ID_nl_strain,SolverID)

% Initialize Global matrices
[K,J]                                                                                  = deal(zeros(n_dofs*n_TotalNodes,n_dofs*n_TotalNodes));                              
[F_internal_Global,F_external_Global,Residual_Global]                                  = deal(zeros(n_dofs*n_TotalNodes,1));              
[storage_strain_loadstep,storage_damage_loadstep,storage_nonlocal_strain_loadstep]     = deal(zeros());

% Partition dofs into displacement and non-local strain
u    = dofs(ID_u);

if SolverID == 2
    e_nl = dofs(ID_nl_strain);
end

% Looping on all elements
for el=1:n_TotalElements
    
    % Identify the Global Coordinates of the nodes of the current element under consideration
    [GlobalCoordinates_CurrentElement] = func_GlobalCoordinates_CurrentElement(n_NodesPerElement,el,connect_LocaltoGlobal_Nodenumbers,Globalcoordinates_EachNode);
    
    % Identify the displacement of the current element under consideration
    u_current_element=[u(el,1);u(el+1,1)];    

    % Identify the non-local strain of the current element under consideration
    if SolverID == 2
        e_nl_current_element = [e_nl(el,1);e_nl(el+1,1)];
    end

    if SolverID == 1
        % Calculate the local stiffness matrix and consistent tangent matrix of the current element under consideration
        [K_el,J_el,f_internal_el,f_external_el,residual_el,storage_damage_loadstep,storage_strain_loadstep] = func_Elstiff_local(ST, ConstitutiveMatrix,weight_integrationpoint,n_NodesPerElement,el,YoungsModulus,n_IntegrationPoints,u_current_element,DamageThreshholdStrain,a,b,GlobalCoordinates_CurrentElement,n_Global,storage_damage_loadstep_con,storage_damage_loadstep,storage_strain_loadstep);
    elseif SolverID == 2
        [K_el,J_el,f_internal_el,f_external_el,residual_el,storage_damage_loadstep,storage_strain_loadstep,storage_nonlocal_strain_loadstep] = func_Elstiff_Gradient(ST, ConstitutiveMatrix,weight_integrationpoint,n_NodesPerElement,el,YoungsModulus,n_IntegrationPoints,u_current_element,e_nl_current_element,DamageThreshholdStrain,a,b,c,GlobalCoordinates_CurrentElement,storage_damage_loadstep_con,storage_strain_loadstep,storage_nonlocal_strain_loadstep,storage_damage_loadstep);
    end

    % Assemble the Global matrices
    [J, K, F_internal_Global, F_external_Global, Residual_Global] = func_GlobalMatrices(n_dofs, n_NodesPerElement,el,connect_LocaltoGlobal_Nodenumbers,J,K,J_el,K_el,f_internal_el,f_external_el,residual_el,F_internal_Global,F_external_Global,Residual_Global);     
  
end

end


