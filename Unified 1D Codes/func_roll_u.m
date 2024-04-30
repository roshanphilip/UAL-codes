function [u] = func_roll_u(n_TotalNodes,u_bar,u_f)

% This function reassembles the Displacement matrix u

% !Update the way this function works (Use flags)!

for i=1:1:n_TotalNodes
    if(i==1)
        u(i,1)=u_bar(1,1);
    elseif i==n_TotalNodes
        u(i,1)=u_bar(2,1);
    else
        u(i,1)=u_f(i-1,1);
    end
end

end