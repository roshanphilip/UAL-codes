function [gamma] = func_calc_gamma(Beta,del_u,del_lamda,delta_u,delta_lamda,q)
% This function calculates the residual of the ArcLength equation g

gamma=((2*delta_u')*(del_u))+(((2*(Beta^2))*delta_lamda*q'*q)*del_lamda);
end