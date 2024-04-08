function [] = func_plot_force_displacement(u_plot,f_plot,trigger,file_name,save_path)
% This function plots a Force-Displacement graph 

if trigger==1
plot(u_plot,f_plot);

% Graph annotation
xlabel('Displacement');
ylabel('Force');

fig_name=file_name(1:end-4);
graph_name='Force vs Displacement';
[t,s] = title(fig_name,graph_name);
t.FontSize = 16;
s.FontAngle = 'italic';

% Save the figure 
image_name=join([fig_name,"_",graph_name]);

% Go to the folder where the images need to be saved
cd (save_path);
% saveas (gcf,image_name);
end
end