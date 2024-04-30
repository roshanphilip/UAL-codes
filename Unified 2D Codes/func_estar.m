function [e_star,s] = func_estar(exx, eyy, gxy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ======================== EQUIVALENT STRAIN E_STAR =======================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the principal strains
e_princ_i = (exx + eyy)/2 + sqrt(((exx-eyy)/2)^2 + (gxy/2)^2);
e_princ_ii = (exx + eyy)/2 - sqrt(((exx-eyy)/2)^2 + (gxy/2)^2);

% Compute the equivalent strain e_star
e_star = sqrt((func_relu(e_princ_i))^2 + (func_relu(e_princ_ii))^2);

%%%%%%%%%%%%%%%% Derivatives %%%%%%%%%%%%%%%%
% The derivatives below are calculated using MATLAB symbolic and the 
% simplify function
if e_star>0
% Derivative of estar wrt exx 
destar_dexx = -(2*((sign(exx/2 + eyy/2 - ((exx/2 - eyy/2)^2 + gxy^2/4)^(1/2))*((exx/2 - eyy/2)/(2*((exx/2 - eyy/2)^2 + gxy^2/4)^(1/2)) - 1/2))/2 + (exx/2 - eyy/2)/(4*((exx/2 - eyy/2)^2 + gxy^2/4)^(1/2)) - 1/4)*(exx/4 + eyy/4 + abs(exx/2 + eyy/2 - ((exx/2 - eyy/2)^2 + gxy^2/4)^(1/2))/2 - ((exx/2 - eyy/2)^2 + gxy^2/4)^(1/2)/2) - 2*((exx/2 - eyy/2)/(4*((exx/2 - eyy/2)^2 + gxy^2/4)^(1/2)) + (sign(exx/2 + eyy/2 + ((exx/2 - eyy/2)^2 + gxy^2/4)^(1/2))*((exx/2 - eyy/2)/(2*((exx/2 - eyy/2)^2 + gxy^2/4)^(1/2)) + 1/2))/2 + 1/4)*(exx/4 + eyy/4 + abs(exx/2 + eyy/2 + ((exx/2 - eyy/2)^2 + gxy^2/4)^(1/2))/2 + ((exx/2 - eyy/2)^2 + gxy^2/4)^(1/2)/2))/(2*((exx/4 + eyy/4 + abs(exx/2 + eyy/2 - ((exx/2 - eyy/2)^2 + gxy^2/4)^(1/2))/2 - ((exx/2 - eyy/2)^2 + gxy^2/4)^(1/2)/2)^2 + (exx/4 + eyy/4 + abs(exx/2 + eyy/2 + ((exx/2 - eyy/2)^2 + gxy^2/4)^(1/2))/2 + ((exx/2 - eyy/2)^2 + gxy^2/4)^(1/2)/2)^2)^(1/2));

% Derivative of estar wrt eyy
destar_deyy = (2*((sign(exx/2 + eyy/2 - ((exx/2 - eyy/2)^2 + gxy^2/4)^(1/2))*((exx/2 - eyy/2)/(2*((exx/2 - eyy/2)^2 + gxy^2/4)^(1/2)) + 1/2))/2 + (exx/2 - eyy/2)/(4*((exx/2 - eyy/2)^2 + gxy^2/4)^(1/2)) + 1/4)*(exx/4 + eyy/4 + abs(exx/2 + eyy/2 - ((exx/2 - eyy/2)^2 + gxy^2/4)^(1/2))/2 - ((exx/2 - eyy/2)^2 + gxy^2/4)^(1/2)/2) - 2*((exx/2 - eyy/2)/(4*((exx/2 - eyy/2)^2 + gxy^2/4)^(1/2)) + (sign(exx/2 + eyy/2 + ((exx/2 - eyy/2)^2 + gxy^2/4)^(1/2))*((exx/2 - eyy/2)/(2*((exx/2 - eyy/2)^2 + gxy^2/4)^(1/2)) - 1/2))/2 - 1/4)*(exx/4 + eyy/4 + abs(exx/2 + eyy/2 + ((exx/2 - eyy/2)^2 + gxy^2/4)^(1/2))/2 + ((exx/2 - eyy/2)^2 + gxy^2/4)^(1/2)/2))/(2*((exx/4 + eyy/4 + abs(exx/2 + eyy/2 - ((exx/2 - eyy/2)^2 + gxy^2/4)^(1/2))/2 - ((exx/2 - eyy/2)^2 + gxy^2/4)^(1/2)/2)^2 + (exx/4 + eyy/4 + abs(exx/2 + eyy/2 + ((exx/2 - eyy/2)^2 + gxy^2/4)^(1/2))/2 + ((exx/2 - eyy/2)^2 + gxy^2/4)^(1/2)/2)^2)^(1/2));
 
% Derivative of estar wrt gxy
destar_dgxy = (2*((gxy*(sign(exx/2 + eyy/2 + (gxy^2 + (exx - eyy)^2)^(1/2)/2) + 1)^2*(exx + eyy + (gxy^2 + (exx - eyy)^2)^(1/2)))/(8*sign(exx/2 + eyy/2 + (gxy^2 + (exx - eyy)^2)^(1/2)/2)*(gxy^2 + (exx - eyy)^2)^(1/2)) - (gxy*(sign(exx/2 + eyy/2 - (gxy^2 + (exx - eyy)^2)^(1/2)/2) + 1)^2*(exx + eyy - (gxy^2 + (exx - eyy)^2)^(1/2)))/(8*sign(exx/2 + eyy/2 - (gxy^2 + (exx - eyy)^2)^(1/2)/2)*(gxy^2 + (exx - eyy)^2)^(1/2))))/((exx + eyy + 2*abs(exx/2 + eyy/2 + (gxy^2 + (exx - eyy)^2)^(1/2)/2) + (gxy^2 + (exx - eyy)^2)^(1/2))^2 + (exx + eyy + 2*abs(exx/2 + eyy/2 - (gxy^2 + (exx - eyy)^2)^(1/2)/2) - (gxy^2 + (exx - eyy)^2)^(1/2))^2)^(1/2);

% Gather the derivatives in vector s
s = [destar_dexx; destar_deyy; destar_dgxy];

if isnan(s)
        s(:) = 0;
end
else
    s = [0.0; 0.0; 0.0];

end
end