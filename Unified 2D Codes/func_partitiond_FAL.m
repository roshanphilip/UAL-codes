function [dofs_E, dofs_F,q,ID_dofs_list_at_ebc,ID_dofs_list_at_nbc,ID_q_dofs,n_essential,n_free] = func_partitiond_FAL(fixnodes_total, dofs,nnodes,Applied_Force_Load,ID_dofs_list_at_ebc,ID_dofs_list_at_nbc,q,trigger)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ===================== PARTITION DISPLACEMENT VECTOR =====================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ndof2 = 2;

if trigger==1
% FAL only
% Partition fixed nodes tensor into dofs with zero applied load and the
% applied force load tensor q

% Identify the location of applied force load nodes
flag_fixed_nodes_FAL=fixnodes_total(3,:)==0;
flag_total_nodes=1:ndof2*nnodes;
ID_q_nodes_FAL_input=find(flag_fixed_nodes_FAL~=1);

% Identify the location of fixed zero nodes
ID_non_q_nodes_FAL=find(flag_fixed_nodes_FAL==1);

fixnodes=fixnodes_total(:,ID_non_q_nodes_FAL);
q_tensor=fixnodes_total(:,ID_q_nodes_FAL_input);

% Identify the essential nodes
for i=1:1:size(ID_non_q_nodes_FAL,2)
    ID_dofs_list_at_ebc(i,1)=fixnodes(1,i)*fixnodes(2,i);
end

% Identify the q nodes
for i=1:1:size(ID_q_nodes_FAL_input,2)
    ID_q_dofs(i,1)=q_tensor(1,i)*q_tensor(2,i);
end

Total_dofs=[1:1:(ndof2*nnodes)]';
ID_dofs_list_at_nbc=setdiff(Total_dofs,ID_dofs_list_at_ebc);


q=zeros(ndof2*nnodes,1);
q(ID_q_dofs,1)=Applied_Force_Load;


% Rearrange the f vector and create the f_E, f_F vectors
% dofs_partitioned = dofs(ID_dofs_list_partitioned); 
dofs_E = dofs(ID_dofs_list_at_ebc);
dofs_F = dofs(ID_dofs_list_at_nbc);

n_essential = size(ID_dofs_list_at_ebc,1);
n_free      = size(ID_dofs_list_at_nbc,1);
else 
    ID_q_dofs=[];
    n_essential=[];
    n_free=[];

    dofs_E = dofs(ID_dofs_list_at_ebc);
    dofs_F = dofs(ID_dofs_list_at_nbc);

end