function [K_el,J_el,f_internal_el,f_external_el,residual_el,storage_damage_loadstep,storage_strain_loadstep] = func_Elstiff_local(ST, ConstitutiveMatrix,weight_integrationpoint,n_NodesPerElement,el,YoungsModulus,n_IntegrationPoints,u_current_element,DamageThreshholdStrain,a,b,GlobalCoordinates_CurrentElement,n_Global,storage_damage_loadstep_con,storage_damage_loadstep,storage_strain_loadstep)

% Initialize local variables
[J_el,K_el]=deal(zeros(n_NodesPerElement,n_NodesPerElement));
[f_internal_el,f_external_el]=deal(zeros(n_NodesPerElement,1));
% Identify the local unit normal vector matrix 
normal_vector=[n_Global(:,el),n_Global(:,el+1)];
Area=1;

% Looping over all the integration points in the current element under consideration
for i=1:1:n_IntegrationPoints

    % Calculate Jacobian of the current element under consideration
    [J] = func_Jacobian(GlobalCoordinates_CurrentElement);

    % Calculate the shape functions of the current element
    [N] = func_ShapeFunctions(n_NodesPerElement);

    % Calculate the derivative of the shape functions of the current element
    [B] = func_ShapeFunctionDerivatives(J);

    % Calculate the equivalent strain
    [Se] = func_EquivalentStrain(B,u_current_element);

    % Calculate damage
    [Damage,dDdestar] = func_damage(DamageThreshholdStrain,a,b,Se);

    % Apply Clausius-Duhem inequality
    if Damage < storage_damage_loadstep_con(1,el) && storage_damage_loadstep_con(1,el) - Damage > ST
        Damage   = storage_damage_loadstep_con(1,el);
        dDdestar = 0;
    end

    % Calculate stress
    [Stress] = func_stress(Damage,Se,YoungsModulus,el);

    % Store the damage in each element during every load cycle
    storage_damage_loadstep(1,el)=Damage;

    % Store the strain in each element during every load cycle
    storage_strain_loadstep(1,el)=Se;

    % Calculate element stiffness matrix 
    for i = 1 : n_NodesPerElement
        for j=1:1:n_NodesPerElement
            K_el(i,j) = K_el(i,j) + ((1-Damage)*B(i)*ConstitutiveMatrix(1,el)*B(j)*weight_integrationpoint(n_IntegrationPoints)*J);
        end
    end

    % Calculate element tangent stiffness matrix
    [J_el] = func_calculate_tangent_stiffness(Damage,dDdestar,B,u_current_element,ConstitutiveMatrix,el,weight_integrationpoint,n_IntegrationPoints,J);

    % Calculate internal force
    f_internal_el = f_internal_el+(B'*Stress*weight_integrationpoint(n_IntegrationPoints)*J);

    % Calculate the reaction force at each element  
    % f_external_el=f_external_el-1*normal_vector*N'*Area*Stress;
    f_external_el = zeros(2,1);

    % Calculate residual vector 
    residual_el = f_internal_el - f_external_el;

end
end

