function [dofs_partitioned, dofs_E, dofs_F, numberof_nodes_at_ebc, numberof_nodes_at_nbc, Nodes_ID_list_at_ebc, Nodes_ID_list_at_nbc] = func_partitiond_fixnodes(fixnodes_applied, dofs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ===================== PARTITION DISPLACEMENT VECTOR =====================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Include global variables
func_include_flags;

% Set up the flags vector. It will indicate the positions of the nodes
% associated with the essential boundary
flags = zeros(ndof*nnodes,1);

% Find the number of nodes at the essential and the natural boundary
numberof_nodes_at_ebc = size(fixnodes_applied,2);
numberof_nodes_at_nbc = ndof*nnodes - numberof_nodes_at_ebc;

% Loop over the number of nodes at the essential boundary and assign a 
% value of 2 in the corresponding position of the [flags] vector
for i = 1:numberof_nodes_at_ebc
    position_of_ebc = ndof*(fixnodes_applied(1,i)-1) + fixnodes_applied(2,i);
    flags(position_of_ebc) = 2;
    dofs(position_of_ebc) = fixnodes_applied(3,i);
end

% Find the nodes IDs at the essential and natural boundary
Nodes_ID_list_at_ebc = find(flags==2);
Nodes_ID_list_at_nbc = find(flags~=2);
Nodes_ID_list_partitioned = [Nodes_ID_list_at_ebc; Nodes_ID_list_at_nbc];

% Rearrange the dofs vector and create the dofs_E, dofs_F vectors
dofs_partitioned = dofs(Nodes_ID_list_partitioned); 
dofs_E = dofs(Nodes_ID_list_at_ebc);
dofs_F = dofs(Nodes_ID_list_at_nbc);

end


