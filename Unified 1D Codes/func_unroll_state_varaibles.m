function [x_u_bar,x_u_f,x_f_reaction_essential,x_m_bar] = func_unroll_state_varaibles(state_variable_vector,M_bar,M)
% This function unrolls the state variable vector into its components 

% Initialize flags matrix
flags=zeros((M+M_bar+M_bar)+1,1);

% Initialize the sub matrices
sub_1=zeros(M_bar,1);
sub_2=zeros(M,1);
sub_3=zeros(M_bar,1);
sub_4=0;

% Assign an identifier to the position of the sub matrices within the flags matrix 
for i=1:1:((M+M_bar+M_bar)+1)
    if i<=M_bar
        flags(i,1)=1;
    elseif i>M_bar && i<=M+M_bar
        flags(i,1)=2;
    elseif i>M+M_bar && i<=(M+M_bar+M_bar)
        flags(i,1)=3;
    else
        flags(i,1)=4;
    end
end

% Store the position of the sub matrices within the flags matrix into another variable
ID_sub_1=find(flags==1);
ID_sub_2=find(flags==2);
ID_sub_3=find(flags==3);
ID_sub_4=find(flags==4);

% Disassemble the converged variable matrix 
sub_1=state_variable_vector(ID_sub_1);
sub_2=state_variable_vector(ID_sub_2);
sub_3=state_variable_vector(ID_sub_3);
sub_4=state_variable_vector(ID_sub_4);

% Assign the values of the submatrices to the corresponding working variables
x_u_bar=sub_1;
x_u_f=sub_2;
x_f_reaction_essential=sub_3;
x_m_bar=sub_4;
end