function [del_f_rct_disp_ebc, del_u_f, del_m_bar,del_e_nl,exit_counter] = func_solve_system_partitioned_UAL_gradient(struct_sub_Matrix_1,R1,R2,R3,R4,M_bar_e,M_e)
% This function calcualtes the corrections based on the partitioned
% consistent scheme

exit_counter = [];

% Calculate del_m_bar
[del_m_bar,A,B,C,D,E,F] = func_calc_del_m_bar_UAL_gradient(struct_sub_Matrix_1{:},R1,R2,R3,R4,M_bar_e,M_e); 

% Calculate del_f_rct
del_f_rct_disp_ebc=E+(F*del_m_bar);

% Calculate del_u_f
del_u_f=C+(D*del_m_bar);

% Calculate del_e_nl
del_e_nl=A+(B*del_m_bar);

end