function [J,K,F_internal_Global,F_external_Global,Residual_Global] = func_GlobalMatrices(n_dofs, n_NodesPerElement,el,connect_LocaltoGlobal_Nodenumbers,J,K,J_el,K_el,f_internal_el,f_external_el,residual_el,F_internal_Global,F_external_Global,Residual_Global);
% This function assembles the Global matrices

% Assemble Consistent Tangent and Stiffness matrix
for a = 1:n_NodesPerElement
    for i = 1:n_dofs
        for b = 1:n_NodesPerElement
            for k = 1:n_dofs
                rw = n_dofs*(connect_LocaltoGlobal_Nodenumbers(a,el)-1)+i;
                cl = n_dofs*(connect_LocaltoGlobal_Nodenumbers(b,el)-1)+k;
                J(rw,cl) = J(rw,cl) + J_el(n_dofs*(a-1)+i,n_dofs*(b-1)+k);
                K(rw,cl) = K(rw,cl) + K_el(n_dofs*(a-1)+i,n_dofs*(b-1)+k);
            end
        end
    end
end

% Assemble Internal force, External force and Residual vector
for a = 1:n_NodesPerElement
    for i = 1:n_dofs
        rw = n_dofs*(connect_LocaltoGlobal_Nodenumbers(a,el)-1)+i;
        F_internal_Global(rw,1) = F_internal_Global(rw,1) + f_internal_el(n_dofs*(a-1)+i,1);
        F_external_Global(rw,1) = F_external_Global(rw,1) + f_external_el(n_dofs*(a-1)+i,1);
        Residual_Global(rw,1)   = Residual_Global(rw,1)   + residual_el(n_dofs*(a-1)+i,1);
    end
end

end

