function [J] = func_Jacobian(GlobalCoordinates_CurrentElement)

    J=(GlobalCoordinates_CurrentElement(1,2)-GlobalCoordinates_CurrentElement(1,1))/2;      % This is J=dx/dT i.e. the mapping between the real and bi-unit domain 

end

