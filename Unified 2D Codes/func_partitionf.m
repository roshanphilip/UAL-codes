function [f_partitioned, f_E, f_F] = func_partitionf(fixnodes_applied, f)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ===================== PARTITION GLOBAL FORCE VECTOR =====================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Include global variables
func_include_flags;

% Construct the force vector based on the essential and free boundaries. 
% [f] is partitioned: f_E, f_F

% Set up the flags vector. It will indicate the positions of the nodes
% associated with the essential boundary
flags = zeros(ndof*nnodes,1);

% Find the number of nodes at the essential and the natural boundary
numberof_nodes_at_ebc = size(fixnodes_applied,2);

% Loop over the number of nodes at the essential boundary and assign a
% value of 2 in the corresponding position of the [flags] vector
for i = 1:numberof_nodes_at_ebc
    position_of_ebc = ndof*(fixnodes_applied(1,i)-1) + fixnodes_applied(2,i);
    flags(position_of_ebc) = 2;
end

% Find the nodes IDs at the essential and natural boundary
Nodes_ID_list_at_ebc = find(flags==2);
Nodes_ID_list_at_nbc = find(flags~=2);
Nodes_ID_list_partitioned = [Nodes_ID_list_at_ebc; Nodes_ID_list_at_nbc];

% Rearrange the f vector and create the f_E, f_F vectors
f_partitioned = f(Nodes_ID_list_partitioned); 
f_E = f(Nodes_ID_list_at_ebc);
f_F = f(Nodes_ID_list_at_nbc);

end
