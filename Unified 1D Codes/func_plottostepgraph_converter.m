function [x2,y2] = func_plottostepgraph_converter(x,y)
%This function takes the input of the material point variable (y-axis) and
%location of gauss points (a-xis) and outputs the coordinates needed to
%plot the step graph 

x=x';
y=y';

% Original coordinates 
size_original=size(x,1);

% Modified coordinates (for step graph)
l_mod_graph=(2*size(x,1))-1;
l_mod_middle=(l_mod_graph+1)/2;

x2=zeros(l_mod_graph,1);
y2=zeros(l_mod_graph,1);

[max_y,pos_max_y]=max(y);
mod_pos_max_y=(pos_max_y*2)-1;

% Calculate the new x positions for the step graph 
j=1;
for i=1:1:size_original
    
    if i~=pos_max_y
        x2(j,1)=x(i,1);
        x2(j+1,1)=x(i,1);
        j=j+2;
    elseif  i==pos_max_y
        x2(j,1)=x(i,1);
        j=j+1;
    end
       
end
clear j 

% Calculate the new y positions for the step graph
[y2] = func_stepfunction_yvalues(y);

x2;
y2;

end

