function [del_f_reaction, del_u_f, del_m_bar,del_e_nl] = func_calc_dels_PC_scheme(Residual_Global, g, delta_u_bar, delta_u_f, delta_f_reaction_essential, J, Beta, Applied_Displacement_Load, ID_prescribed_u, ID_free_u, f_rct)
% This function finds the values of corrections for f_reaction, u_f and m_bar

% Partitional Residual_Global vector
r_e   = Residual_Global(ID_prescribed_u,1) + f_rct;
r_f   = Residual_Global(ID_free_u,1);

% Partition consistent tangent matrix
Jff   = J(ID_free_u,ID_free_u);
Jfe   = J(ID_free_u,ID_prescribed_u);
Jef   = J(ID_prescribed_u,ID_free_u);
Jee   = J(ID_prescribed_u,ID_prescribed_u);

% Calculate g1,g2,g3
g1    = 2 * delta_u_bar';    
g2    = 2 * delta_u_f';    
g3    = 2 * (Beta^2) * delta_f_reaction_essential';

% Calculate del_m_bar
del_u_NR   = -Jff\r_f;
del_u_g    = -Jff\(Jfe * Applied_Displacement_Load);
num        = g + (g2 * del_u_NR) - (g3 * ((Jef * del_u_NR) + r_e));
den        = (g1 * Applied_Displacement_Load) + (g2 * del_u_g) - (g3 * ((Jee * Applied_Displacement_Load) + (Jef * del_u_g)));

del_m_bar  = ((-1)*(num))/(den);

% Calculate del_u_f
del_u_f    = del_u_NR + (del_m_bar*del_u_g);

% Calculate del_f_rct
del_f_reaction = (-1) * (r_e+(del_m_bar * Jee * Applied_Displacement_Load) + (Jef * del_u_f));

del_e_nl   = zeros;

end