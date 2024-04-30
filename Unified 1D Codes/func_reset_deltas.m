function [delta_u_bar,delta_u_f,delta_f_reaction_essential,delta_m_bar] = func_reset_deltas(convergance_flag,iteration_flag,current_deltas,converged_deltas,M,M_bar,Applied_Displacement_Load)
% This function resets the values of all the deltas based on the current state of the system  
    
% Unroll converged deltas
[delta_u_bar_con,delta_u_f_con,delta_f_reaction_essential_con,delta_m_bar_con] = func_unroll_state_varaibles(converged_deltas,M_bar,M);

% Unroll current deltas 
[delta_u_bar,delta_u_f,delta_f_reaction_essential,delta_m_bar] = func_unroll_state_varaibles(current_deltas,M_bar,M);


% A. If converged in the previous loadstep, set the deltas to the last converged values of deltas from the previous loadstep  
if convergance_flag==1 && iteration_flag==0
    delta_m_bar=delta_m_bar_con;
    delta_u_f=delta_u_f_con;
    delta_u_bar=delta_u_bar_con;
    delta_f_reaction_essential=delta_f_reaction_essential_con;
    
% B. If convergence not reached in the last loadstep after k_max iterations, reduce the last converged deltas by a factor of x (obtained from eq 37)
elseif convergance_flag==0 && iteration_flag==0
    delta_m_bar=delta_m_bar_con;
    delta_u_bar=delta_m_bar*Applied_Displacement_Load;
    delta_u_f=delta_u_f_con;
    delta_f_reaction_essential=delta_f_reaction_essential_con;

% C. If convergence not reached in the last loadstep but k<k_max, keep the deltas as they are 
elseif convergance_flag==0 && iteration_flag==1
    delta_m_bar;
    delta_u_f;
    delta_u_bar;
    delta_f_reaction_essential;
end

end