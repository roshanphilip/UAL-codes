function [] = func_plot_strain(strain_solution,n_TotalElements,l_Total,trigger,file_name,save_path)
% This function plots the Damage across the domain

if trigger==1
% Calcualte the Gauss point locations %
l_EachElement=l_Total/n_TotalElements;
nodal_locations=zeros(1,n_TotalElements);
nodal_locations(1,1)=l_EachElement/2;
for i=2:1:n_TotalElements
    nodal_locations(1,i)=nodal_locations(1,i-1)+l_EachElement;
end

% Assign x-axis values 
x_0=nodal_locations;

% Assign y-axis values
y_0=strain_solution;

% Calculate the step graph coordinates 
[x,y] = func_plottostepgraph_converter(x_0,y_0);

% Plot the damage across the domain
plot(x,y);

%----------Graph annotations----------%
% X-axis
xlabel ('Gauss point location');
xlim([x(1,1) x(size(x,1),1)]);

% Y-axis
ylabel ('Strain');

% Figure name and style
fig_name=file_name(1:end-4);
graph_name='Strain across the domain';
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