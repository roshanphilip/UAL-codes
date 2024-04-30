function [] = func_plot(u_plot,f_plot,storage_damage_loadstep_global_matrix,storage_strain_loadstep_global_matrix,storage_nonlocal_strain_loadstep_global_matrix,n_TotalElements,l_Total,plot_flag,save_flag,file_path,save_path,SolverID,file_name)
% This function saves the parameters needed for the graphs and plots the graphs if need be

plot_storage(:,1) = u_plot;
plot_storage(:,2) = f_plot;

damage_solution=storage_damage_loadstep_global_matrix(size(storage_damage_loadstep_global_matrix,1),:);
strain_solution=storage_strain_loadstep_global_matrix(size(storage_strain_loadstep_global_matrix,1),:);
nonlocal_strain_solution = storage_nonlocal_strain_loadstep_global_matrix(size(storage_damage_loadstep_global_matrix,1),:);

% Choose which graphs to plot  % 0 - Don't plot graph, 1 - Plot graph
if SolverID == 1
    graph_flag=[1,1,1,0];            % [Force-Displacement, Damage, Strain,Non-local strain]
elseif SolverID == 2
    graph_flag=[1,1,1,0];            % [Force-Displacement, Damage, Strain,Non-local strain]
end

% Save key variables from the solution into a new file 
if save_flag==1
    cd (save_path)
    save (file_name)
end

% Plot graphs 
% Return to the current directory 
cd (file_path)
if plot_flag==1
    for i=1:size(graph_flag,2)
        if i==1 && graph_flag(1,i)==1
            % Plot the Force-Displacement graph
            figure(i)
            cd (file_path)
            func_plot_force_displacement(u_plot,f_plot,graph_flag(1,i),file_name,save_path);
            
        elseif i==2 && graph_flag(1,i)==1
            % Plot the Damage across the domain
            figure(i)
            cd (file_path)
            func_plot_damage(damage_solution,n_TotalElements,l_Total,graph_flag(1,i),file_name,save_path);

        elseif i==3 && graph_flag(1,i)==1
            % Plot the Strain across the domain
            figure(i)
            cd (file_path)
            func_plot_strain(strain_solution,n_TotalElements,l_Total,graph_flag(1,i),file_name,save_path);
        elseif i==4 && graph_flag(1,i)==1
            % Plot the Strain across the domain
            figure(i)
            cd (file_path)
            func_plot_strain(nonlocal_strain_solution,n_TotalElements,l_Total,graph_flag(1,i),file_name,save_path);
        end
    end
end
end
