function [x,y] = func_separate_x_y_UAL(dofs)
% This function separates the x and y components of the input variable
% Assuming that the input variable is in the format:
% Input = [x1;y1;x2;y2;........xn;yn]
n=1;
m=1;
size_dofs=size(dofs,1);
for i=1:1:size_dofs
    if rem(i,2)~=0
        x(n,1)=dofs(i,1);
        n=n+1;
    else 
        y(m,1)=dofs(i,1);
        m=m+1;
    end
end

end