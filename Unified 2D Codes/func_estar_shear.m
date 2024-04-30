function [e_star,s] = func_estar(exx, eyy, gxy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ======================== EQUIVALENT STRAIN E_STAR =======================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the principal strains
e_princ_i = (exx + eyy)/2 + sqrt(((exx-eyy)/2)^2 + (gxy/2)^2);
e_princ_ii = (exx + eyy)/2 - sqrt(((exx-eyy)/2)^2 + (gxy/2)^2);

% Compute the equivalent strain e_star
e_star = (0.75*exx) + (0.75*eyy) + (17787670415741550*exx^2 + 29711231904181648*exx*eyy + 17787670415741550*eyy^2 + 1466027231825363*gxy^2)^(1/2)/167772160;
 

%%%%%%%%%%%%%%%% Derivatives %%%%%%%%%%%%%%%%
% The derivatives below are calculated using MATLAB symbolic and the 
% simplify function
if e_star>0
% Derivative of estar wrt exx 
destar_dexx = (8893835207870775*exx + 7427807976045412*eyy + 62914560*(17787670415741550*exx^2 + 29711231904181648*exx*eyy + 17787670415741550*eyy^2 + 1466027231825363*gxy^2)^(1/2))/(83886080*(17787670415741550*exx^2 + 29711231904181648*exx*eyy + 17787670415741550*eyy^2 + 1466027231825363*gxy^2)^(1/2));

% Derivative of estar wrt eyy
destar_deyy = (7427807976045412*exx + 8893835207870775*eyy + 62914560*(17787670415741550*exx^2 + 29711231904181648*exx*eyy + 17787670415741550*eyy^2 + 1466027231825363*gxy^2)^(1/2))/(83886080*(17787670415741550*exx^2 + 29711231904181648*exx*eyy + 17787670415741550*eyy^2 + 1466027231825363*gxy^2)^(1/2));

% Derivative of estar wrt gxy
destar_dgxy = (1466027231825363*gxy)/(167772160*(17787670415741550*exx^2 + 29711231904181648*exx*eyy + 17787670415741550*eyy^2 + 1466027231825363*gxy^2)^(1/2));


% Gather the derivatives in vector s
s = [destar_dexx; destar_deyy; destar_dgxy];

if isnan(s)
    s(:) = 0;
end
else
    s = [0.0; 0.0; 0.0];
end
end