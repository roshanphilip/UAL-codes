function [g] = func_UAL_arclength_eq(delta_u_bar,delta_u_f,delta_f_reaction_essential,delta_e_nl,ArcLength,Beta,SolverID)
% This function computes residual of the Arc length equation

if SolverID == 1     % Local damage
    g1_i   =  delta_u_bar'*delta_u_bar;
    g2_i   =  delta_u_f'*delta_u_f;
    g3_i   =  (Beta^2)*(delta_f_reaction_essential'*delta_f_reaction_essential);
    g5_i   =  -((ArcLength)^2);
    g      =  g1_i + g2_i + g3_i + g5_i;

elseif SolverID == 2 % Gradient damage
    g1_i   =  delta_u_bar'*delta_u_bar;
    g2_i   =  delta_u_f'*delta_u_f;
    g3_i   =  (Beta^2)*(delta_f_reaction_essential'*delta_f_reaction_essential);
    g4_i   =  (delta_e_nl' * delta_e_nl);
    g5_i   =  -((ArcLength)^2);
    g      =  g1_i+g2_i+g3_i+g4_i+g5_i;
end

end