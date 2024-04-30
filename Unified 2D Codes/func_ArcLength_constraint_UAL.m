function [g] = func_ArcLength_constraint_UAL(SolverID, delta_e_nl, delta_u_bar,delta_u_f,delta_f_reaction_essential,ArcLength,Beta)
% This function computes the Arc length equation
% For a spherical contraint equation (Using eq 18.)

if SolverID == 1 % Local damage
    % Initialize g,g1,g2,g3
    [g,g1_i,g2_i,g3_i]=deal(0);

    g1_i = delta_u_bar' * delta_u_bar;
    g2_i = delta_u_f'   * delta_u_f;
    g3_i = (Beta^2)     * (delta_f_reaction_essential' * delta_f_reaction_essential);
    g4_i = -((ArcLength)^2);

    g    = g1_i + g2_i + g3_i + g4_i;

elseif SolverID == 2 % Non-local damage
    % Initialize g,g1,g2,g3
    [g,g1_i,g2_i,g3_i]=deal(0);
    
    g1_i = delta_u_bar' * delta_u_bar;
    g2_i = delta_u_f'   * delta_u_f;
    g3_i = delta_e_nl'  * delta_e_nl;
    g4_i = (Beta^2)     * (delta_f_reaction_essential' * delta_f_reaction_essential);
    g5_i = -((ArcLength)^2);

    g    = g1_i + g2_i + g3_i + g4_i + g5_i;

end


end