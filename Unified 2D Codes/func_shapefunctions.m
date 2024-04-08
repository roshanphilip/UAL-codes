function [N] = func_shapefunctions(xi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ============================ SHAPE FUNCTIONS ============================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculates shape functions for various element types
% Specify just the special case: 2D linear quadrilateral elements
N = zeros(4,1);
N(1) = 0.25*(1.-xi(1))*(1.-xi(2));
N(2) = 0.25*(1.+xi(1))*(1.-xi(2));
N(3) = 0.25*(1.+xi(1))*(1.+xi(2));
N(4) = 0.25*(1.-xi(1))*(1.+xi(2));


end