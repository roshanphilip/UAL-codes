function [x] = func_cosine_limits_UAL(x)
% This function ensures that the variable has a value between -1 and 1.
% acosd gives a complex output if the values are out of this range

if x>1
    x=1;
elseif x<-1
    x=-1;
else 
    x=x;
end

end