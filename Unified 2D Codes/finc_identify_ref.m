function [ref, ref_node_ID, ref_node_pos] = finc_identify_ref(fixnodes, ndof, model_name)
% This function identfies the reference dof to measure the displacement applied at each increment in the case of NR and FAL routines 

% Identify the direction of load application 
if model_name == 'SNS'
    direction_load = 1;    % In the Single Notch Shear problem, load is applied in the postive x direction 
elseif model_name=='SSNT_Coarse' || model_name=='SSNT_Fine' || model_name=='TNT_Coarse' || model_name=='TNT_Fine' || model_name=='SNT' 
    direction_load = 2;
end

% Identify the first node with a displacement in the +ve Y direction
% (assuming displacement is applied along the Y axis)
for i=1:1:size(fixnodes,2)
    if fixnodes(3,i)>0 && fixnodes(2,i) == direction_load 
        ref_node_ID  = fixnodes(1,i);
        ref_node_pos = i; 
        break
    end
end

ref = (ref_node_ID - 1) * ndof + direction_load;

end