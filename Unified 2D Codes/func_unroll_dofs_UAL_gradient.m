function [u_bar,u_f,e_star] = func_unroll_dofs_UAL_gradient(dofs,ID_dofs_list_u_bar,ID_dofs_list_u_f,ID_dofs_list_nl_strain)
% This function unrolls dofs 

u_f=dofs(ID_dofs_list_u_f,1);
u_bar=dofs(ID_dofs_list_u_bar,1);
e_star=dofs(ID_dofs_list_nl_strain,1);

end