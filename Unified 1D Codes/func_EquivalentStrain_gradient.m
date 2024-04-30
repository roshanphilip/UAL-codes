function [Se,S,dEeq_dEz,Se_non_local] = func_EquivalentStrain_gradient(B,u_current_element,N,e_nl_current_element)
% This function calculates the equivalent strain in each element
                     
S=0;
Se=0;

% Flip the elements in the u_current_element so that the strain get calculated correctly
% If this is not done, we will always get a -ve value for strain which doesn't make phyisical sense for the problem we are trying to solve here 

S=B*u_current_element;                      % Strain in the current element 
Se=Se+S;                                    % Equivalent strain 

dEeq_dEz=1;                                 % dEeq/dEz

Se_non_local=N*e_nl_current_element;

end

