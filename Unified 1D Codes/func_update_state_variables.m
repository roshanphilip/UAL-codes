function [u_bar,u_f,m_bar,u,f_reaction_essential,current_state_variables] = func_update_state_variables(M_bar,M,converged_variables,current_deltas,n_TotalNodes,Applied_Displacement_Load)
% This function updates the state variables based on the last converged
% values + delta of the current iteration 

% Unroll current deltas
[delta_u_bar,delta_u_f,delta_f_reaction_essential,delta_m_bar] = func_unroll_state_varaibles(current_deltas,M_bar,M);

% Unroll converged state variables
[u_bar_con,u_f_con,f_reaction_essential_con,m_bar_con] = func_unroll_state_varaibles(converged_variables,M_bar,M);

% Update the state variables
m_bar=m_bar_con+delta_m_bar;
u_bar=m_bar*Applied_Displacement_Load;
u_f=u_f_con+delta_u_f;
u = func_roll_u(n_TotalNodes,u_bar,u_f);
f_reaction_essential=f_reaction_essential_con+delta_f_reaction_essential;

% Roll current state variable values 
[current_state_variables] = func_roll_state_variables(u_bar,u_f,f_reaction_essential,m_bar,M_bar,M);

end