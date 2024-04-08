function [state_variable_vector] = func_roll_state_variables(x_u_bar,x_u_f,x_f_reaction_essential,x_m_bar,M_bar,M)
% This function assembles the state variables into a vector

% Initialize flags matrix
flags=zeros((M+M_bar+M_bar)+1,1);

% Initialize the converged variables matrix
state_variable_vector=zeros((M+M_bar+M_bar)+1,1);

% A. Assembling the converged_variables matrix

% Identify the sub matrices
sub_1=x_u_bar;
sub_2=x_u_f;
sub_3=x_f_reaction_essential;
sub_4=x_m_bar;

% Assign an identifier to the position of the sub matrices within the flags matrix
for i=1:1:((2*(M_bar))+1+M)
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

% Assemble the converged variable matrix
state_variable_vector(ID_sub_1)=sub_1;
state_variable_vector(ID_sub_2)=sub_2;
state_variable_vector(ID_sub_3)=sub_3;
state_variable_vector(ID_sub_4)=sub_4;

end
