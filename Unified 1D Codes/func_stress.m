function [Stress] = func_stress(Damage,Se,YoungsModulus,el)

% This function calculates the stress at each gauss point 

Stress=(1-Damage)*YoungsModulus(el)*Se;

end

