function [Bu, Be] = func_Bmatrix(dNdx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ===================== STRAIN-DISPLACEMENT MATRIX (B) ====================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% NON LOCAL GRADIENT IMPLEMENTATION %%%%% 
% Defines the strain-displacement matrix B for each element point, by appropriately expanding the dNdx input argument. 
% Bu (3x8) is the strain-displacement matrix used in the calculation of ux, uy      
% Be (2x4) is the strain-displacement matrix used in the calculation of enonlocal

% Specify just the special case: 2D linear quadrilateral elements. 
% Transpose the dNdx matrix (for 2D quadrilateral convert from 4x2 to 2x4)
dNdx_transp = dNdx';

Bu = zeros(3,8);
Bu(1,1) = dNdx_transp(1,1);
Bu(2,2) = dNdx_transp(2,1);
Bu(3,1) = dNdx_transp(2,1);
Bu(3,2) = dNdx_transp(1,1);
Bu(1,3) = dNdx_transp(1,2);
Bu(2,4) = dNdx_transp(2,2);
Bu(3,3) = dNdx_transp(2,2);
Bu(3,4) = dNdx_transp(1,2);
Bu(1,5) = dNdx_transp(1,3);
Bu(2,6) = dNdx_transp(2,3);
Bu(3,5) = dNdx_transp(2,3);
Bu(3,6) = dNdx_transp(1,3);
Bu(1,7) = dNdx_transp(1,4);
Bu(2,8) = dNdx_transp(2,4);
Bu(3,7) = dNdx_transp(2,4);
Bu(3,8) = dNdx_transp(1,4);

Be = dNdx_transp;

end






