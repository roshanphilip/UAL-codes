function [r_e,r_f,g] = func_sort_residuals(residual,M,M_bar)
% This function calculates the residuals at the current position of the
% system

% Initialize variables
r_e=zeros(M_bar,1);
r_f=zeros(M,1);
g=0;

for i=1:1:M_bar+M+1
    if i>=1 && i<=M_bar
        r_e(i,1)=residual(i,1);
    elseif i>=M_bar && i<=M_bar+M
        r_f(i-M_bar,1)=residual(i,1);
    elseif i==M_bar+M+1
        g=residual(i,1);
    end
end