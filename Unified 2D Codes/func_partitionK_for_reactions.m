function [K_E, K_EF, K_FE, K_F, numberof_nodes_at_ebc, numberof_nodes_at_nbc, Nodes_ID_list_at_ebc, Nodes_ID_list_at_nbc] = func_partitionK_for_reactions(fixnodes_applied, K)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =================== PARTITION GLOBAL STIFFNESS MATRIX ===================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Include global variables
func_include_flags;

% Construct the global stiffness matrix based on the essential and free 
% boundaries. [K] is partitioned: K_E, K_EF, K_FE, K_F 

% Set up the flags vector. It will indicate the positions of the nodes
% associated with the essential boundary
flags = zeros(ndof2*nnodes,1);

% Find the number of nodes at the essential and the natural boundary
numberof_nodes_at_ebc = size(fixnodes_applied,2);
numberof_nodes_at_nbc = ndof2*nnodes - numberof_nodes_at_ebc;

% Loop over the number of nodes at the essential boundary and assign a
% value of 2 in the corresponding position of the [flags] vector
for i = 1:numberof_nodes_at_ebc
    position_of_ebc = ndof2*(fixnodes_applied(1,i)-1) + fixnodes_applied(2,i);
    flags(position_of_ebc) = 2;
end

% Find the nodes IDs at the essential and natural boundary
Nodes_ID_list_at_ebc = find(flags==2);
Nodes_ID_list_at_nbc = find(flags~=2);

% Rearrange the K matrix and create the K_E, K_EF, K_FE, K_F matrices
%K_partitioned = K(Nodes_ID_list_partitioned,Nodes_ID_list_partitioned); 
K_E  = K(Nodes_ID_list_at_ebc, Nodes_ID_list_at_ebc);
K_EF = K(Nodes_ID_list_at_ebc, Nodes_ID_list_at_nbc);
K_FE = K(Nodes_ID_list_at_nbc, Nodes_ID_list_at_ebc);
K_F  = K(Nodes_ID_list_at_nbc, Nodes_ID_list_at_nbc);

% % [fixnodes] contains prescribed displacements both at the essential
% % bounday (3rd row value is zero) and potentially at the natural boundary
% % (3rd row value is non-zero, for displacement-control load).
% fixnodes_E = fixnodes(:,(fixnodes(3,:)==0));



end