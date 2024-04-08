function [del_m_bar,A,B,C,D,E,F] = func_calc_del_m_bar_UAL_gradient(sub_1,sub_2,sub_3,sub_4,sub_5,sub_6,sub_7,sub_8,sub_9,sub_10,sub_11,sub_12,sub_13,sub_14,sub_15,sub_16,R1,R2,R3,R4,M_bar_e,M_e)
% This function calculates the parameters used in the partitioned
% consistent scheme

% f_rct parameters 
del_f_rct_1=-R1;
del_f_rct_2=-sub_2;
del_f_rct_3=-sub_3;
del_f_rct_4=-sub_4;

% u_f parameters
inv_sub_6=inv(sub_6);
del_u_f_1=(-1)*(inv_sub_6)*R2;
del_u_f_2=(-1)*(inv_sub_6)*(sub_7);
del_u_f_3=(-1)*(inv_sub_6)*(sub_8);

% m_bar parameters
del_m_bar_1=(-1)*(R3/sub_11);
del_m_bar_2=(-1)*(sub_9/sub_11);
del_m_bar_3=(-1)*(sub_10/sub_11);

% e_nl parameters
inv_sub_16=inv(sub_16);
del_e_nl_1=(-1)*(inv_sub_16)*R4;
del_e_nl_2=(-1)*(inv_sub_16)*(sub_14);
del_e_nl_3=(-1)*(inv_sub_16)*(sub_15);

% del_e_nl parameters
A_den=(eye(M_bar_e+M_e,M_bar_e+M_e)-(del_e_nl_2*del_u_f_3));
A=A_den\(del_e_nl_1+(del_e_nl_2*del_u_f_1));
B=A_den\(del_e_nl_3+(del_e_nl_2*del_u_f_2));

% del_u_f parameters
C=del_u_f_1+(del_u_f_3*A);
D=del_u_f_2+(del_u_f_3*B);

% del_f_rct parameters
E=del_f_rct_1+(del_f_rct_2*C)+(del_f_rct_4*A);
F=(del_f_rct_2*D)+del_f_rct_3+(del_f_rct_4*B);

% Calculate del_m
del_m_bar=(del_m_bar_1+(del_m_bar_2*E)+(del_m_bar_3*C))/(1-(del_m_bar_2*F)-(del_m_bar_3*D));

end