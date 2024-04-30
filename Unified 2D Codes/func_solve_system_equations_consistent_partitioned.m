function [del_f_rct_disp_ebc, del_u_f, del_m_bar,exit_counter] = func_solve_system_equations_consistent_partitioned(J_F,J_E,J_EF,J_FE,Res_F_E_rct,Res_F_F,g_constraint,delta_u_bar,delta_u_f,delta_f_rct_disp_ebc,listofnodes_ebc,listofnodes_nbc,Beta,Applied_Displacement_Load)
% This function finds the values of corrections for f_reaction, u_f and m_bar

exit_counter=[];

% Initialize del_dofs
del_dofs=zeros((size(listofnodes_ebc,1)+size(listofnodes_nbc,1)),1);

% Calculate g1,g2,g3
g1=2*delta_u_bar';    
g2=2*delta_u_f';    
g3=2*(Beta^2)*delta_f_rct_disp_ebc';

% Calculate del_m_bar
del_u_NR=(-1)*(J_F\Res_F_F);
del_u_g=(-1)*(J_F\(J_FE*Applied_Displacement_Load));
num=g_constraint+(g2*del_u_NR)-(g3*((J_EF*del_u_NR)+Res_F_E_rct));
den=(g1*Applied_Displacement_Load)+(g2*del_u_g)-(g3*((J_E*Applied_Displacement_Load)+(J_EF*del_u_g)));

del_m_bar=((-1)*(num))/(den);

del_u_bar=del_m_bar*Applied_Displacement_Load;

% Calculate del_u_f
del_u_f=del_u_NR+(del_m_bar*del_u_g);

% Calculate del_f_rct
del_f_rct_disp_ebc=(-1)*(Res_F_E_rct+(del_m_bar*J_E*Applied_Displacement_Load)+(J_EF*del_u_f));

% Assemble del_dofs
del_dofs(listofnodes_ebc) = del_u_bar;
del_dofs(listofnodes_nbc) = del_u_f;

end