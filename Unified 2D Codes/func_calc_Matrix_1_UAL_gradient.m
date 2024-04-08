function [Matrix_1,struct_sub_Matrix_1] = func_calc_Matrix_1_UAL_gradient(Jff,Jpp,Jpf,Jfp,Kue_p,Kue_f,Keu_p,Keu_f,Kee,Applied_Displacement_Load,Beta,delta_u_bar,delta_u_f,delta_f_reaction_essential,M_bar_u,M_u,M_bar_e,M_e)
% This function calculates Matrix_1

% Initialize Matrix_1
Matrix_1=zeros((M_bar_u+M_u+M_bar_e+M_e+1),(M_bar_u+M_u+M_bar_e+M_e+1));

sub_1=eye(M_bar_u,M_bar_u);
sub_2=Jpf;
sub_3=Jpp*Applied_Displacement_Load;
sub_4=Kue_p;
sub_5=zeros(M_u,M_bar_u);
sub_6=Jff;
sub_7=Jfp*Applied_Displacement_Load;
sub_8=Kue_f;
sub_9=2*(Beta^2)*delta_f_reaction_essential';
sub_10=2*delta_u_f';
sub_11=2*delta_u_bar'*Applied_Displacement_Load;
sub_12=zeros(1,M_bar_e+M_e);
sub_13=zeros(M_bar_e+M_e,M_bar_u);
sub_14=Keu_f;
sub_15=Keu_p*Applied_Displacement_Load;
sub_16=Kee;

[Matrix_1] = func_assemble_Matrix_1_UAL_gradient(M_bar_u,M_u,M_bar_e,M_e,sub_1,sub_2,sub_3,sub_4,sub_5,sub_6,sub_7,sub_8,sub_9,sub_10,sub_11,sub_12,sub_13,sub_14,sub_15,sub_16);
struct_sub_Matrix_1={sub_1,sub_2,sub_3,sub_4,sub_5,sub_6,sub_7,sub_8,sub_9,sub_10,sub_11,sub_12,sub_13,sub_14,sub_15,sub_16};

end