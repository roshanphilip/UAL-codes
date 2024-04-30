function [xi] = func_integrationpoints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ===================== INTEGRATION POINTS POSITIONS ======================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify just the special case: 2D linear quadrilateral elements
xi = zeros(2,4);
xi(1,1) = -0.5773502692;
xi(2,1) = xi(1,1);
xi(1,2) = -xi(1,1);
xi(2,2) = xi(1,1);
xi(1,3) = xi(1,1);
xi(2,3) = -xi(1,1);
xi(1,4) = -xi(1,1);
xi(2,4) = -xi(1,1);

end