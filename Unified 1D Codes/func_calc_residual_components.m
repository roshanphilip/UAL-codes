function [g,Jff,Jee,Jef,Jfe,F_f,F_e,storage_damage_loadstep,storage_strain_loadstep] = func_calc_residual_components(struct_0,current_deltas,Beta,M_bar,M,ArcLength)
% This function calculates the components of the residual matrix and the variables needed to calculate the dels in the next iteration 

% Unroll the current deltas 
[delta_u_bar,delta_u_f,delta_f_reaction_essential,delta_m_bar] = func_unroll_state_varaibles(current_deltas,M_bar,M);

% Calculate g
[g] = func_constraint_equation_18(delta_u_bar,delta_u_f,delta_f_reaction_essential,ArcLength,Beta,1);

% Calculate f_internal, Tangent stiffness components and damage across the domain 
[~,~,~,~,Jff,Jee,Jef,Jfe,F_e,F_f,~,~,storage_damage_loadstep,storage_strain_loadstep] = func_CalculateGlobalMatrices(struct_0{:});


end