function [residual_norm,residual] = func_calc_residual_norm(f_reaction_essential,F_int_essential,F_int_free,g,M,M_bar)
% This function calculates the residuals at the current position of the
% system

% Calculate residual matrix 
residual=zeros(M_bar+M+1,1);
r_e=f_reaction_essential+F_int_essential;

for i=1:1:M_bar+M+1
    if i>=1 && i<=M_bar
        residual(i,1)=r_e(i,1);
    elseif i>=M_bar && i<=M_bar+M
        residual(i,1)=F_int_free(i-M_bar,1);
    elseif i==M_bar+M+1
        residual(i,1)=g;
    end
end

% Calculate norm of the residual 
residual_norm=norm(residual,2);

% Calculate residuals of each of the state variables 
clc
res_norm_f_rct=norm(r_e,2);
res_norm_u_f=norm(F_int_free,2);
res_norm_g=norm(g,2);

end
