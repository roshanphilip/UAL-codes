function [convergance_flag,iteration_flag,converged_variables,converged_deltas,ArcLength,n,storage_damage_loadstep_global_matrix,storage_strain_loadstep_global_matrix,storage_norm_loadstep_global_matrix,u_plot,f_plot,p] = func_check_convergence_ArcLength_local(current_state_variables,current_deltas,residual_norm,tolerance,k,k_max,k_des,ArcLength,ArcLength_0,Gamma,n,storage_damage_loadstep,storage_strain_loadstep,storage_damage_loadstep_global_matrix,storage_strain_loadstep_global_matrix,storage_norm_loadstep_global_matrix,converged_variables,converged_deltas,u_plot,f_plot,M_bar,M,p,Applied_Displacement_Load)
% This function checks convergence and either moves to the next iteration or next loadstep

% A. Convergence is reached 
if residual_norm<=tolerance 
    
    % Unroll converged variables
    [u_bar,u_f,f_reaction_essential,m_bar] = func_unroll_state_varaibles(current_state_variables,M_bar,M);

    % Store the value of the converged variables and converged deltas 
    converged_variables=current_state_variables;
    converged_deltas=current_deltas;

    % Update ArcLength 
    if k<5
        if u_bar>Applied_Displacement_Load*0.8
            ArcLength=min(1e-4,(10^(log10(ArcLength)+0.2)));
        else 
            ArcLength=min(1e-4,(10^(log10(ArcLength)+0.2)));
        end
    elseif k>=5 && k<=12
        ArcLength=ArcLength;
    elseif k>12
        ArcLength=max(1e-24,(10^(log10(ArcLength)-0.2)));
    end

    % Save the current damage and strain across the domain 
    storage_damage_loadstep_global_matrix(n,:)=storage_damage_loadstep;
    storage_strain_loadstep_global_matrix(n,:)=storage_strain_loadstep;
    storage_norm_loadstep_global_matrix(n,k)=residual_norm;
    
    % Update increment counter
    n=n+1;

    % Update iteration and convergence flag
    convergance_flag=1;
    iteration_flag=0;

    % Save the F-displacement parameters
    u_plot(p,1)=u_bar(2,1);
    f_plot(p,1)=f_reaction_essential(1,1);
    p=p+1;
    
% B. Convergence isn't reached but k<k_max (Go to the next iteration)
elseif residual_norm>tolerance && k<k_max

    % Update iteration and convergence flag
    convergance_flag=0;
    iteration_flag=1;

% C. Convergence isn't reached and k=k_max (We need to go back to the last converged position of the system) 
elseif residual_norm>tolerance && k==k_max

   % Update ArcLength 
    if k<5
        if u_bar>Applied_Displacement_Load*0.8
            ArcLength=min(1e-4,(10^(log10(ArcLength)+0.2)));
        else 
            ArcLength=min(1e-4,(10^(log10(ArcLength)+0.2)));
        end
    elseif k>=5 && k<=12
        ArcLength=ArcLength;
    elseif k>12
        ArcLength=max(1e-24,(10^(log10(ArcLength)-0.2)));
    end
    
    % Update iteration and convergence flag
    convergance_flag=0;
    iteration_flag=0;

end