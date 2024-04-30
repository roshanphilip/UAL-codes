function [Jpf,Jpp,Jff,Jfp,Jue_p,Jue_f,Jeu_p,Jeu_f,Jee] = func_sort_stiffness_UAL_gradient(J,ID_dofs_list_u_p,ID_dofs_list_u_f,ID_dofs_list_nl_strain,ID_dofs_list_disp,fixnodes_applied)
% This function partitions th estiffness matrix

% This function sorts the global matrices into the submatrices 

% Split the J matrix based on the displacements and non local strains
Juu=J(ID_dofs_list_disp,ID_dofs_list_disp);
Jue=J(ID_dofs_list_disp,ID_dofs_list_nl_strain);
Jeu=J(ID_dofs_list_nl_strain,ID_dofs_list_disp);
Jee=J(ID_dofs_list_nl_strain,ID_dofs_list_nl_strain);


% Find the ID of the free and essential dofs
flags=zeros(size(J,2),1);
flags_strain=zeros(size(Jee,2),1);
flags(ID_dofs_list_u_p,1)=1;
dofs_ID_free=find(flags~=1);
dofs_ID_all=find(flags>-1);
dofs_ID_strain_all=find(flags_strain==0);

% Sort J matrix based on prescribed and free displacements
Jpp=J(ID_dofs_list_u_p,ID_dofs_list_u_p);
Jpf=J(ID_dofs_list_u_p,ID_dofs_list_u_f);
Jfp=J(ID_dofs_list_u_f,ID_dofs_list_u_p);
Jff=J(ID_dofs_list_u_f,ID_dofs_list_u_f);

% Find the position of the prescribed and free nodes within the Jue and Jeu
% matrix 
flags_Ke=zeros(size(Jue,1),1);
for i = 1:size(fixnodes_applied,2)
    position_of_ebc_e = 2*(fixnodes_applied(1,i)-1) + fixnodes_applied(2,i);
    flags_Ke(position_of_ebc_e) = 2;
end
ID_dofs_list_u_bar_e=find(flags_Ke==2);
ID_dofs_list_u_f_e=find(flags_Ke~=2);

% Sort Jue matrix
Jue_p=Jue(ID_dofs_list_u_bar_e,dofs_ID_strain_all);
Jue_f=Jue(ID_dofs_list_u_f_e,dofs_ID_strain_all);

% Sort Jeu matrix
Jeu_p=Jeu(dofs_ID_strain_all,ID_dofs_list_u_bar_e);
Jeu_f=Jeu(dofs_ID_strain_all,ID_dofs_list_u_f_e);

end