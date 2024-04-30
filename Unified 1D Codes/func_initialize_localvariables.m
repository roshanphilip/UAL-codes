function [K_local,K_local_elastic,J_analytical_local,J11,J12,J21,J22,f_internal_u_el,f_rct_el,normal_vector] = func_initialize_localvariables(n_NodesPerElement,n_Global,el)
% Initialize the value of K_local and J_analytical_local (and its componenets) to zero for each element

K_local=zeros(n_NodesPerElement,n_NodesPerElement);
K_local_elastic=zeros(n_NodesPerElement,n_NodesPerElement);
J_analytical_local=zeros(n_NodesPerElement,n_NodesPerElement);
[f_internal_u_el,f_rct_el]=deal(zeros(n_NodesPerElement,1));
[J11,J12,J21,J22]=deal(0);

% Identify the local unit normal vector matrix
normal_vector=[n_Global(:,el),n_Global(:,el+1)];
end

