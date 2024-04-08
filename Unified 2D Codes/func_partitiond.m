function [dofs_partitioned, dofs_E, dofs_F,ID_dofs_list_partitioned,ID_dofs_list_at_ebc,ID_dofs_list_at_nbc] = func_partitiond(fixnodes_applied, dofs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ===================== PARTITION DISPLACEMENT VECTOR =====================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Include global variables
func_include_flags;

% Set up the flags vector. It will indicate the positions of the nodes
% associated with the essential boundary
flags = zeros(ndof*nnodes,1);

% Find the number of nodes at the essential boundary
numberof_nodes_at_ebc = size(fixnodes_applied,2);

% Loop over the number of nodes at the essential boundary and assign a
% value of 2 in the corresponding position of the [flags] vector
for i = 1:numberof_nodes_at_ebc
    position_of_ebc = ndof*(fixnodes_applied(1,i)-1) + fixnodes_applied(2,i);
    flags(position_of_ebc) = 2;
end

% Find the nodes IDs at the essential and natural boundary
ID_dofs_list_at_ebc = find(flags==2);
ID_dofs_list_at_nbc = find(flags~=2);
ID_dofs_list_partitioned = [ID_dofs_list_at_ebc; ID_dofs_list_at_nbc];

% Rearrange the f vector and create the f_E, f_F vectors
dofs_partitioned = dofs(ID_dofs_list_partitioned); 
dofs_E = dofs(ID_dofs_list_at_ebc);
dofs_F = dofs(ID_dofs_list_at_nbc);


end