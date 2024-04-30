function [struct_sub_Matrix_1] = func_calc_Matrix_1(J,Applied_Displacement_Load,Beta,delta_u_bar,delta_u_f,delta_f_reaction_essential,M_bar,M,ID_prescribed_u,ID_free_u,ID_nl_strain)
% This function calculates Matrix_1

% Partition J
Jee        =   J(ID_nl_strain,ID_nl_strain);
Jeu_p      =   J(ID_nl_strain,ID_prescribed_u);
Jue_p      =   J(ID_prescribed_u,ID_nl_strain);
Jue_f      =   J(ID_free_u,ID_nl_strain);
Jeu_f      =   J(ID_nl_strain,ID_free_u);

Jff        =   J(ID_free_u,ID_free_u);
Jfp        =   J(ID_free_u,ID_prescribed_u);
Jpf        =   J(ID_prescribed_u,ID_free_u);
Jpp        =   J(ID_prescribed_u,ID_prescribed_u);

% Matrix_1 submatrices
sub_1      =   eye(M_bar,M_bar);
sub_2      =   Jpf;
sub_3      =   Jpp*Applied_Displacement_Load;
sub_4      =   Jue_p;
sub_5      =   zeros(M,M_bar);
sub_6      =   Jff;
sub_7      =   Jfp*Applied_Displacement_Load;
sub_8      =   Jue_f;
sub_9      =   2*(Beta^2)*delta_f_reaction_essential';
sub_10     =   2*delta_u_f';
sub_11     =   2*delta_u_bar'*Applied_Displacement_Load;
sub_12     =   zeros(1,M_bar+M);
sub_13     =   zeros(M_bar+M,M_bar);
sub_14     =   Jeu_f;
sub_15     =   Jeu_p*Applied_Displacement_Load;
sub_16     =   Jee;

struct_sub_Matrix_1={sub_1,sub_2,sub_3,sub_4,sub_5,sub_6,sub_7,sub_8,sub_9,sub_10,sub_11,sub_12,sub_13,sub_14,sub_15,sub_16};

end