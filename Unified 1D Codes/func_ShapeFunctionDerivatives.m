function [B] = func_ShapeFunctionDerivatives(J)

% // This function calculates the derivatives of the shape functions of the current element under consideration //
    dNdT=[-0.5,0.5];   % This is dN/dT which is the derivative of the shape function w.r.t. the coordinates of the bi-unit domain 
    
    dNdx=dNdT*(1/J);   % This is (dN/dT)*(dT/dx)=(dN/dT)*(1/J)=dN/dx
    
    B=dNdx;

end

