function [residual_norm,R_f] = func_calc_residual(F_f,F_e,n_TotalNodes)
% This function calculates the norm of the residual

% Calculate residuals
R_e=zeros(2,1);
R_f=F_f;

% Roll residual vector
[Global_Residual] = func_roll_u(n_TotalNodes,R_e,R_f);

% Calculate residual norm 
residual_norm=norm(Global_Residual,2);

end