function [dNdxi] = func_shapefunctionderivs(xi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====================== SHAPE FUNCTION DERIVATIVES =======================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify just the special case: 2D linear quadrilateral elements
dNdxi      = zeros(4,2);
dNdxi(1,1) = -0.25*(1.-xi(2));
dNdxi(1,2) = -0.25*(1.-xi(1));
dNdxi(2,1) = 0.25*(1.-xi(2));
dNdxi(2,2) = -0.25*(1.+xi(1));
dNdxi(3,1) = 0.25*(1.+xi(2));
dNdxi(3,2) = 0.25*(1.+xi(1));
dNdxi(4,1) = -0.25*(1.+xi(2));
dNdxi(4,2) = 0.25*(1.-xi(1));

end