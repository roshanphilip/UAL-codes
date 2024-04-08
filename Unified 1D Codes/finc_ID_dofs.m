function [ID_prescribed_u,ID_free_u,ID_u,ID_nl_strain,ID_free, q] = finc_ID_dofs(n_dofs, n_TotalNodes, Applied_Force_Load, dofs, SolverID, RoutineID)
% This function identifies the ID of dofs based on different
% classifications

% For UAL and NR
if RoutineID == 1 || RoutineID == 3   
    q = [];
    if SolverID == 1   % Local damage
        ID_prescribed_u = [1 ; size(dofs,1)];
        ID_u            = (1 : size(dofs,1))';
        ID_free_u       = setdiff(ID_u,ID_prescribed_u);
        ID_free         = ID_free_u;
        ID_nl_strain    = zeros; 
    elseif SolverID == 2 % Gradient damage
        ID_prescribed_u = [1 ; (size(dofs,1)-1)];
        flags           = zeros(size(dofs,1),1);
        for i = 2 : 2 : size(dofs,1)
            flags(i,1)  = 3;
        end
        ID_nl_strain    = find(flags == 3);
        ID_u            = find(flags ~= 3);
        ID_free_u       = setdiff(ID_u,ID_prescribed_u);
        ID_free         = union(ID_nl_strain,ID_free_u);
    end
% For FAL
elseif RoutineID == 2
    ID_prescribed_u     = [1];
    if SolverID == 1   % Local damage
        ID_u            = (1 : size(dofs,1))';
        ID_free_u       = setdiff(ID_u,ID_prescribed_u);
        ID_free         = ID_free_u;
        ID_nl_strain    = zeros;
    elseif SolverID == 2 % Gradient damage
        flags           = zeros(size(dofs,1),1);
        for i = 2 : 2 : size(dofs,1)
            flags(i,1)  = 3;
        end
        ID_nl_strain    = find(flags == 3);
        ID_u            = find(flags ~= 3);
        ID_free_u       = setdiff(ID_u,ID_prescribed_u);
        ID_free         = union(ID_nl_strain,ID_free_u);
    end
    % Calculate applied force vector q 
    q                          = zeros(n_dofs * n_TotalNodes,1);
    q(ID_free_u(end,1),1)      = Applied_Force_Load;     % Applied force vector (FAL only)
end