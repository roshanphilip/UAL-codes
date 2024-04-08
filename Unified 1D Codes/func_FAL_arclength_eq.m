function [gamma] = func_FAL_arclength_eq(Beta,del_dofs,del_lamda,delta_dofs,delta_lamda,q,ArcLength)
% This function calculates the residual of the ArcLength equation Gamma

gamma=(delta_dofs+del_dofs)' * (delta_dofs+del_dofs) + (Beta^2) * (2 * (delta_lamda + del_lamda)) * (q' *q) - ArcLength^2;

end