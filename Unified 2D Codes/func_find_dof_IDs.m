function [ID_dofs_list_u_p,ID_dofs_list_u_f,ID_dofs_list_nl_strain,ID_dofs_list_disp,ID_dofs_list_x,ID_dofs_list_y,ID_free_nodes_e,ID_prescribed_nodes_e] = func_find_dof_IDs(fixnodes,ndof,dofs_stored,SolverID)
% This function identify the ID of the u_bar,u_f and non-local equivalent
% strain from the dofs vector

if SolverID == 1 % For local model
    [ID_dofs_list_nl_strain,ID_dofs_list_disp,ID_dofs_list_x,ID_dofs_list_y,ID_free_nodes_e,ID_prescribed_nodes_e]=deal([]);
    
    for i = 1:size(fixnodes,2)
        position_of_ebc = ndof*(fixnodes(1,i)-1) + fixnodes(2,i);
        flags(position_of_ebc) = 2;
    end

    % Identify ID of dofs
    ID_dofs_list_u_f=find(flags==0);
    ID_dofs_list_u_p=find(flags==2);

elseif SolverID == 2 % For non-local gradient  model
    flags=zeros(size(dofs_stored,1),1);

    for i = 1:size(fixnodes,2)
        position_of_ebc = ndof*(fixnodes(1,i)-1) + fixnodes(2,i);
        flags(position_of_ebc) = 2;
    end

    for i=3:3:size(dofs_stored,1)
        flags(i,1)=3;
    end

    % Identify ID of dofs
    ID_dofs_list_nl_strain = find(flags==3);
    ID_dofs_list_u_f=find(flags==0);
    ID_dofs_list_u_p=find(flags==2);
    ID_dofs_list_disp=find(flags~=3);

    % Identify position of x and y displacement dofs
    flag_xy=zeros(size(dofs_stored,1),1);
    x_i=0;
    y_i=1;

    for i=1:3:size(dofs_stored,1)
        position_of_x=x_i+i;
        flag_xy(position_of_x)=1;
        position_of_y=y_i+i;
        flag_xy(position_of_y)=2;
    end

    % Identify ID of dofs
    ID_dofs_list_x=find(flag_xy==1);
    ID_dofs_list_y=find(flag_xy==2);

    % ID of the free and prescribed nodes (for nl_strain calculations only)
    ID_prescribed_nodes_e=unique(fixnodes(1,:));
    X=[1:size(ID_dofs_list_nl_strain,1)];
    ID_free_nodes_e=X;
    ID_free_nodes_e(ID_prescribed_nodes_e)=[];

    ID_prescribed_nodes_e=ID_prescribed_nodes_e';
    ID_free_nodes_e=ID_free_nodes_e';


end