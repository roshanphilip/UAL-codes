function [Kff,Kee,Kef,Kfe,Jff,Jee,Jef,Jfe,F_e,F_f,u_e,u_f] = func_SortingGlobalMatrices(M_u_p,M_u_f,pos_essential_dofs,J_tangent,K,F_internal_Global,F_reaction_Global,u)
% This function sorts the Global matrices 
% 1. Sort J_tangent 

% Initialize flag vector 
flags=zeros(M_u_p+M_u_f,1);                            

% Assign a value of 1 to the position of essential and free dofs in the flags vector
flags(pos_essential_dofs)=1;

% Find the ID of the free and essential dofs
dofs_ID_essential=find(flags==1);
dofs_ID_free=find(flags~=1);

% Sort J_tangent matrix 
Jee=J_tangent(dofs_ID_essential,dofs_ID_essential);
Jef=J_tangent(dofs_ID_essential,dofs_ID_free);
Jfe=J_tangent(dofs_ID_free,dofs_ID_essential);
Jff=J_tangent(dofs_ID_free,dofs_ID_free);

% Sort K matrix 
Kee=K(dofs_ID_essential,dofs_ID_essential);
Kef=K(dofs_ID_essential,dofs_ID_free);
Kfe=K(dofs_ID_free,dofs_ID_essential);
Kff=K(dofs_ID_free,dofs_ID_free);

% 2. Sort F_internal_Global
F_int_e=F_internal_Global(dofs_ID_essential,1);
F_int_f=F_internal_Global(dofs_ID_free,1);

% 3. Sort F_reaction_Global
F_rct_e=F_reaction_Global(dofs_ID_essential,1);
F_rct_f=F_reaction_Global(dofs_ID_free,1);

% Calculate residual forces on the free and essential nodes
F_e=F_int_e+F_rct_e;
F_f=F_int_f;

% 4. Sort displacement 
u_e=u(dofs_ID_essential,1);
u_f=u(dofs_ID_free,1);
                                 
end

