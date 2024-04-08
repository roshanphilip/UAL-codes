function [Se] = func_EquivalentStrain(B,u_current_element)

% This function calculates the equivalent strain in each element
                     
S=0;
Se=0;

S=B*u_current_element;              % Strain in the current element 
Se=Se+S;                                    % Equivalent strain 

end

