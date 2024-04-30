function [u_f,e_star] = func_unroll_free_dofs_UAL_gradient(ID_dofs_list_u_p,ID_dofs_list_u_f,ID_dofs_list_nl_strain,free_dofs,prescribed_dofs)
% This function unrolls the free dofs

dofs_temp=zeros((size(ID_dofs_list_u_p,1)+size(ID_dofs_list_u_f,1)+size(ID_dofs_list_nl_strain,1)),1);

merged_ID=[ID_dofs_list_u_f;ID_dofs_list_nl_strain];
ID_dofs_list_free_dofs=sort(merged_ID);

dofs_temp(ID_dofs_list_u_p,1)=prescribed_dofs;
dofs_temp(ID_dofs_list_free_dofs,1)=free_dofs;

[u_bar,u_f,e_star] = func_unroll_dofs_UAL_gradient(dofs_temp,ID_dofs_list_u_p,ID_dofs_list_u_f,ID_dofs_list_nl_strain);

end