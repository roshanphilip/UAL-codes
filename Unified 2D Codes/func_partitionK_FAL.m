function [K_E, K_EF, K_FE, K_F] = func_partitionK_FAL(J, ID_dofs_list_at_ebc, ID_dofs_list_at_nbc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =================== PARTITION GLOBAL STIFFNESS MATRIX ===================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Include global variables
func_include_flags;

% Rearrange the J matrix and create the K_E, K_EF, K_FE, K_F matrices
%K_partitioned = J(Nodes_ID_list_partitioned,Nodes_ID_list_partitioned); 
K_E  = J(ID_dofs_list_at_ebc, ID_dofs_list_at_ebc);
K_EF = J(ID_dofs_list_at_ebc, ID_dofs_list_at_nbc);
K_FE = J(ID_dofs_list_at_nbc, ID_dofs_list_at_ebc);
K_F  = J(ID_dofs_list_at_nbc, ID_dofs_list_at_nbc);

end