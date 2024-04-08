function [n,p,u_con,u_plot,f_plot,storage_damage_loadstep_global_matrix,storage_strain_loadstep_global_matrix,storage_norm_loadstep_global_matrix] = func_check_convergence(k,k_max,n,p,u,F_e,n_TotalNodes,tolerance,u_con,u_plot,f_plot,storage_damage_loadstep,storage_strain_loadstep,residual_norm,storage_damage_loadstep_global_matrix,storage_strain_loadstep_global_matrix,storage_norm_loadstep_global_matrix)
% This function checks the onvergence at the end of a loadstep 

% Unroll u
[u_e,~]=func_unroll_u(u,n_TotalNodes);
      
% Store the value of the converged variables and converged deltas
u_con=u;

% Save the current damage and strain across the domain
storage_damage_loadstep_global_matrix(n,:)=storage_damage_loadstep;
storage_strain_loadstep_global_matrix(n,:)=storage_strain_loadstep;
storage_norm_loadstep_global_matrix(n,k)=residual_norm;

% Update increment counter
n=n+1;

% Save the F-displacement parameters
u_plot(p,1)=u_e(2,1);
f_plot(p,1)=F_e(2,1);
p=p+1;

