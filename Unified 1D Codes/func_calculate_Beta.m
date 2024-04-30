function [Beta] = func_calculate_Beta(ST, constraint_type,M_bar,n_dofs, c, n_TotalElements,n_TotalNodes,n_NodesPerElement,connect_LocaltoGlobal_Nodenumbers,Globalcoordinates_EachNode,n_IntegrationPoints,ConstitutiveMatrix,weight_integrationpoint,DamageThreshholdStrain,a,b,n_Global,storage_damage_loadstep_con,YoungsModulus,dofs,ID_u,ID_nl_strain,ID_prescribed_u,SolverID)
% This function calculates the value of Beta based on the type of
% constraint chosen 

% A. Shperical constraint (Eq. 19)
if constraint_type==1
    % Calculate consistent tangent matrix J
    [~,J,~,~,~,~,~,~] = func_CalculateGlobalMatrices(ST, n_dofs, c, n_TotalElements,n_TotalNodes,n_NodesPerElement,connect_LocaltoGlobal_Nodenumbers,Globalcoordinates_EachNode,n_IntegrationPoints,ConstitutiveMatrix,weight_integrationpoint,DamageThreshholdStrain,a,b,n_Global,storage_damage_loadstep_con,YoungsModulus,dofs,ID_u,ID_nl_strain,SolverID);
    Jpp = J(ID_prescribed_u,ID_prescribed_u);

    Beta=inv((1/M_bar)*trace(Jpp));

% B. Cylindrical constraint 
elseif constraint_type==2
    Beta=0;
 
% C. Fixed Beta 
elseif constraint_type==3
    Beta=10e-8;
end

end