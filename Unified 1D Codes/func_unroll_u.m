function [u_e,u_f]=func_unroll_u(u,n_TotalNodes)

% This function seperates the free and essential node components from the
% displacement matrix

u_f=zeros(n_TotalNodes-2,1);

for i=1:1:n_TotalNodes    
    if i==1 
        u_e(1,1)=u(i,1);
    elseif i==n_TotalNodes
        u_e(2,1)=u(i,1);
    else 
        u_f(i-1,1)=u(i,1);
    end       
end

end

