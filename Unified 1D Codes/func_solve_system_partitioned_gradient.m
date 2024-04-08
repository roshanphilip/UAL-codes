function [del_f_reaction, del_u_f, del_m_bar,del_e_nl] = func_solve_system_partitioned_gradient(struct_sub_Matrix_1,Residual_Global, g, M_bar, M, ID_prescribed_u, ID_free_u, ID_u, ID_nl_strain, ID_free, f_rct)
% This function calculates the corrections based on the partitioned
% consistent scheme

R1   = Residual_Global(ID_prescribed_u,1) + f_rct;
R2   = Residual_Global(ID_free_u,1);
R3   = g; 
R4   = Residual_Global(ID_nl_strain,1);

% Calculate del_m_bar
[del_m_bar,A,B,C,D,E,F] = func_calc_del_m_bar(struct_sub_Matrix_1{:},R1,R2,R3,R4,M_bar,M); 

% Calculate del_f_rct
del_f_reaction=E+(F*del_m_bar);

% Calculate del_u_f
del_u_f=C+(D*del_m_bar);

% Calculate del_e_nl
del_e_nl=A+(B*del_m_bar);

end