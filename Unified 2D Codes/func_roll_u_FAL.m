function [u] = func_roll_u_FAL(pos_essential_dofs,u_bar,u_f)
% This function reassembles the Displacement matrix u

% Initialize flag vector
flags=zeros(((size(u_bar,1))+(size(u_f,1))),1);

% Assign a value of 1 to the position of essential and free dofs in the flags vector
flags(pos_essential_dofs)=1;

% Initialize the converged variables matrix
u=zeros(((size(u_bar,1))+(size(u_f,1))),1);

% Find the ID of the free and essential dofs
dofs_ID_essential=find(flags==1);
dofs_ID_free=find(flags~=1);

% Assemble the displacement vector 
u(dofs_ID_essential,1)=u_bar;
u(dofs_ID_free,1)=u_f;
end