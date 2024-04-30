function [K_el,J_el,f_internal_el,f_external_el,residual_el,storage_damage_loadstep,storage_strain_loadstep,storage_nonlocal_strain_loadstep] = func_Elstiff_Gradient(ST, ConstitutiveMatrix,weight_integrationpoint,n_NodesPerElement,el,YoungsModulus,n_IntegrationPoints,u_current_element,e_nl_current_element,DamageThreshholdStrain,a,b,c,GlobalCoordinates_CurrentElement,storage_damage_loadstep_con,storage_strain_loadstep,storage_nonlocal_strain_loadstep,storage_damage_loadstep)
 
% This is the Element stiffness matrix for the Displacement controlled ArcLength method with Gradient damage
%============================================================================================%
    % Initialize local variables
    [f_internal_u_el,f_internal_e_nl_el] = deal(zeros(2,1));
    [Juu_el,Jue_el,Jeu_el,Jee_el]        = deal(zeros(2,2));
 
    % Looping over all the integration points in the current element under consideration
    for i=1:1:n_IntegrationPoints

%============================================================================================%
    % Calculating the parameters needed to calculate K_local and J_Analytical_local

    % Calculate Jacobian of the current element under consideration
    [J] = func_Jacobian(GlobalCoordinates_CurrentElement);
    % Calculate the shape functions of the current element
    [N] = func_ShapeFunctions(n_NodesPerElement);
    % Calculate the derivative of the shape functions of the current element
    [B] = func_ShapeFunctionDerivatives(J);
    % Calculate the equivalent strain
    [Se,S,dEeq_dEz,Se_non_local] = func_EquivalentStrain_gradient(B,u_current_element,N,e_nl_current_element);
    
    % Calculate damage
    [Damage,dDdestar] = func_damage(DamageThreshholdStrain,a,b,Se_non_local);
    
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
    % Store the non-local strain in each element during every load cycle 
    storage_nonlocal_strain_loadstep(1,el)=Se_non_local;
        
%============================================================================================%
    % Calculate Tangent stiffness matrix (Gradient)          
    Juu_el    = Juu_el + (B' * (1-Damage) * ConstitutiveMatrix(1,el) * B * weight_integrationpoint(n_IntegrationPoints) * J);   
    Jue_el    = Jue_el + ((-1) * B' * ConstitutiveMatrix(1,el) * dDdestar * S * N * weight_integrationpoint(n_IntegrationPoints) * J);    
    Jeu_el    = Jeu_el + ((-1) * N' * dEeq_dEz * B * weight_integrationpoint(n_IntegrationPoints) * J);   
    Jee_el    = Jee_el + (((N' * N) + (B' * c * B)) * weight_integrationpoint(n_IntegrationPoints) * J);
%============================================================================================%      
    % Assemble consistent tangent matrix (J_el)
    J_el     = [Juu_el Jue_el; Jeu_el Jee_el];

    % Identify elastic stiffness matrix
    K_el     = [Juu_el zeros(2,2); zeros(2,2) zeros(2,2)];
%============================================================================================%  
    % Calculate the internal force (displacement) at each element 
    f_internal_u_el      = f_internal_u_el + B'*Stress*weight_integrationpoint(n_IntegrationPoints)*J;          

    % Calculate the internal force (non-local strain) at each element 
    f_internal_e_nl_el   = f_internal_e_nl_el + ((N'*N*e_nl_current_element) + (B'*c*B*e_nl_current_element) - (N'*Se))*weight_integrationpoint(n_IntegrationPoints)*J;

    % Assemble internal force vector 
    f_internal_el        = [f_internal_u_el(1,1);f_internal_e_nl_el(1,1);f_internal_u_el(2,1);f_internal_e_nl_el(2,1)];
    % External force vector 
    f_external_el        = zeros(4,1);

    % Calculate residual 
    residual_el          = f_internal_el - f_external_el;  

%============================================================================================%  
    % Rearrange consistent tangent and residual 

    % Set up the indices order
    % [u1_1, u2_1, u1_2, u2_2] --> [u1_1, u1_2, u2_1, u2_2]
    index_series = [1 3 2 4];

    % Partition the matrix and vectors to follow the global system 
    J_el          = J_el(index_series,index_series);
    K_el          = K_el(index_series,index_series);
        
    end          
end

