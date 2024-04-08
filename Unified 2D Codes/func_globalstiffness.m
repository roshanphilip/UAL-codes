function [K, J, DRdu, Res_F, F_int, F_ext, history_var_mat, gausspoints_prop_mat, nodes_prop_mat, strain_var_mat] = func_globalstiffness(model_name,dofs,Delastic,history_var_mat_previousinc,num_elem_at_node,n_hood,weights,strain_tolerance,strain_mat_previousinc,IsProj,RoutineID)     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====================== ASSEMBLE THE TANGENT MATRIX ======================
% ===================== ASSEMBLE THE RESIDUAL VECTOR ======================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Include global variables
func_include_flags;

% Variables strain_var_mat and strain_mat_previousinc are only active in
% the UAL local damage case 
if RoutineID == 2 || RoutineID == 3 || SolverID == 2
    strain_var_mat = [];
    strain_mat_previousinc = zeros(size(history_var_mat_previousinc,1),size(history_var_mat_previousinc,2));
end

% Calculate the local damage matrix for the nonlocal integral method
if SolverID == 3
    [damage_mat_it] = func_localdamagemat(dofs); 
end

% -------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SOLVER ROUTINE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
if IsProj == 0

    % Create empty entries for element/nodal properties 
    gausspoints_prop_mat = [];
    nodes_prop_mat = [];
   
    % Initializing global stiffness, global internal force and gather matrices 
    % at zero:
    J                   = zeros(ndof*nnodes,ndof*nnodes);
    K                   = zeros(ndof2*nnodes,ndof2*nnodes);
    DRdu                = zeros(ndof*nnodes,ndof*nnodes);
    Res_F               = zeros(ndof*nnodes,1);
    F_int               = zeros(ndof*nnodes,1);
    F_ext               = zeros(ndof*nnodes,1);
    residual_elpos      = zeros(maxnodes*ndof,maxnodes*ndof);
    residual_elneg      = zeros(maxnodes*ndof,maxnodes*ndof);
    lmncoord            = zeros(ncoord,maxnodes);
    lmndof              = zeros(ndof,maxnodes);
    history_var_mat     = zeros(nelem,4);
    
    % Loop over all the elements
    for lmn = 1:nelem
    
        % Extract coords of nodes, DOF for the current element
        for a = 1:nelnodes(lmn)
            for i = 1:ncoord
                lmncoord(i,a) = coords(i,connect(a,lmn));
            end
            for i = 1:ndof
                lmndof(i,a) = dofs(ndof*(connect(a,lmn)-1)+i);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ANALYTICAL CALCULATION - TANGENT MATRIX (for each element): 
        
        if TangentID == 1

            % -------------------------------------------------------------
            if SolverID == 1 
                [k_el, j_el, history_var_mat(lmn,:), strain_var_mat(lmn,:),~, ~, Res_F_el, f_internal_el, f_external_el] = func_elstif_Local(model_name,lmncoord,lmndof,Delastic,history_var_mat_previousinc(lmn,:),alpha_val,beta_val,e_delta,dmax,n_hood,weights,strain_tolerance,strain_mat_previousinc(lmn,:),1,IsProj,RoutineID);                
            % -------------------------------------------------------------
            elseif SolverID == 2
                [k_el, j_el, history_var_mat(lmn,:), ~, ~, Res_F_el, f_internal_el, f_external_el] = func_elstif_Nonlocgradient(model_name,lmncoord,lmndof,Delastic,history_var_mat_previousinc(lmn,:),g,alpha_val,beta_val,e_delta,dmax,strain_tolerance,n_hood,weights,1,IsProj);
            % -------------------------------------------------------------
            elseif SolverID == 3
                [j_el, history_var_mat(lmn,:), ~, ~, Res_F_el, f_external_el] = func_elstif_Nonlocintegral(lmn,lmncoord,lmndof,Delastic,history_var_mat_previousinc(lmn,:),alpha_val,beta_val,e_delta,dmax,n_hood,weights,damage_mat_it,1,IsProj);
            else
                disp("Check your SolverID - globalstiffness")
            end

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % NUMERICAL APPROXIMATION - TANGENT MATRIX (for each element):
    
        if TangentID == 2
                
            % -------------------------------------------------------------
            if SolverID == 1
                [~, ~, history_var_mat(lmn,:), strain_var_mat(lmn,:), ~, ~, Res_F_el, f_internal_el, f_external_el] = func_elstif_Local(model_name,lmncoord,lmndof,Delastic,history_var_mat_previousinc(lmn,:),alpha_val,beta_val,e_delta,dmax,n_hood,weights,strain_tolerance,strain_mat_previousinc(lmn,:),1,IsProj,RoutineID);    
            % -------------------------------------------------------------
            elseif SolverID == 2
                [~, ~, history_var_mat(lmn,:), ~, ~, Res_F_el, f_internal_el, f_external_el] = func_elstif_Nonlocgradient(model_name,lmncoord,lmndof,Delastic,history_var_mat_previousinc(lmn,:),g,alpha_val,beta_val,e_delta,dmax,strain_tolerance,n_hood,weights,1,IsProj);
            % -------------------------------------------------------------
            elseif SolverID == 3
                [~, history_var_mat(lmn,:), ~, ~, Res_F_el, f_external_el] = func_elstif_Nonlocintegral(lmn,lmncoord,lmndof,Delastic,history_var_mat_previousinc(lmn,:),alpha_val,beta_val,e_delta,dmax,n_hood,weights,damage_mat_it,1,IsProj);
            % -------------------------------------------------------------
            else 
                disp("Check your SolverID - globalstiffness")
            end

            % element stiffness (based on residual at plus/minus lmndof - central difference)
            dl = 1e-10; % Tolerance 
            col_id = 0; % Counter of columns in dRdu (range: 1-12)
        
            for a = 1:nelnodes(lmn)
                for i = 1:ndof
                    
                    % Increase counter
                    col_id = col_id + 1;
        
                    % Compute residual_el at (u+dl(dof)) of size: 12x1, and store 
                    % it in the next column of residual
                    lmndofi = lmndof;    
                    lmndofi(i,a) = lmndofi(i,a) + dl;
                    % -----------------------------------------------------
                    if SolverID == 1
                        [~, ~, ~, ~, ~, ~, residual_elpos(:,col_id), f_internal_el, f_external_el] = func_elstif_Local(model_name,lmncoord,lmndofi,Delastic,history_var_mat_previousinc(lmn,:),alpha_val,beta_val,e_delta,dmax,n_hood,weights,strain_tolerance,strain_mat_previousinc(lmn,:),0,IsProj,RoutineID);
                    % -----------------------------------------------------
                    elseif SolverID == 2
                        [~, ~, ~, ~, ~, residual_elpos(:,col_id), f_internal_el, f_external_el] = func_elstif_Nonlocgradient(model_name,lmncoord,lmndofi,Delastic,history_var_mat_previousinc(lmn,:),g,alpha_val,beta_val,e_delta,dmax,strain_tolerance,n_hood,weights,0,IsProj); % Final size: 12x12
                    % -----------------------------------------------------
                    elseif SolverID == 3
                        [~, ~, ~, ~, residual_elpos(:,col_id), f_external_el] = func_elstif_Nonlocintegral(lmn,lmncoord,lmndofi,Delastic,history_var_mat_previousinc(lmn,:),alpha_val,beta_val,e_delta,dmax,n_hood,weights,damage_mat_it,0,IsProj);
                    % -----------------------------------------------------
                    else
                        disp("Check your SolverID - globalstiffness")
                    end

                    % Compute residual_el at (u-dl(dof)) of size: 12x1, and store 
                    % it in the next column of residual
                    lmndofi = lmndof;    
                    lmndofi(i,a) = lmndofi(i,a) - dl;
                    % -----------------------------------------------------
                    if SolverID == 1 
                        [~, ~, ~, ~, ~, ~, residual_elpos(:,col_id), f_internal_el, f_external_el] = func_elstif_Local(model_name,lmncoord,lmndofi,Delastic,history_var_mat_previousinc(lmn,:),alpha_val,beta_val,e_delta,dmax,n_hood,weights,strain_tolerance,strain_mat_previousinc(lmn,:),0,IsProj,RoutineID);
                    % -----------------------------------------------------
                    elseif SolverID == 2
                        [~, ~, ~, ~, ~, residual_elpos(:,col_id), f_internal_el, f_external_el] = func_elstif_Nonlocgradient(model_name,lmncoord,lmndofi,Delastic,history_var_mat_previousinc(lmn,:),g,alpha_val,beta_val,e_delta,dmax,strain_tolerance,n_hood,weights,0,IsProj); % Final size: 12x12
                    % -----------------------------------------------------
                    elseif SolverID == 3
                        [~, ~, ~, ~, residual_elpos(:,col_id), f_external_el] = func_elstif_Nonlocintegral(lmn,lmncoord,lmndofi,Delastic,history_var_mat_previousinc(lmn,:),alpha_val,beta_val,e_delta,dmax,n_hood,weights,damage_mat_it,0,IsProj);
                    % -----------------------------------------------------
                    else
                        disp("Check your SolverID - globalstiffness")
                    end
                end
            end
 
            % Compute "partial R over partial u" for all dofs at one step 
            dRdu_el = 1/(2*dl) * (residual_elpos - residual_elneg); % Size: 12x12

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ASSEMBLE GLOBAL MATRICES: 
        % a) analytical tangent K matrix
        % b) numerical tangent DRdu matrix
        % c) residual vector
        
        % a) Analytical tangent matrix
        if TangentID == 1
            for a = 1:nelnodes(lmn)
                for i = 1:ndof
                    for b = 1:nelnodes(lmn)
                        for k = 1:ndof
                            rw = ndof*(connect(a,lmn)-1)+i;
                            cl = ndof*(connect(b,lmn)-1)+k;
                            J(rw,cl) = J(rw,cl) + j_el(ndof*(a-1)+i,ndof*(b-1)+k);
                        end
                    end
                end
            end
        end
    
        % b) Numerical tangent matrix
        if TangentID  == 2
            for a = 1:nelnodes(lmn)
                for i = 1:ndof
                    for b = 1:nelnodes(lmn)
                        for k = 1:ndof
                            rw = ndof*(connect(a,lmn)-1)+i;
                            cl = ndof*(connect(b,lmn)-1)+k;
                            DRdu(rw,cl) = DRdu(rw,cl) + dRdu_el(ndof*(a-1)+i,ndof*(b-1)+k);
                        end
                    end
                end
            end
        end

         % c) Stiffness matrix 
        if TangentID == 1
            for a = 1:nelnodes(lmn)
                for i = 1:ndof2
                    for b = 1:nelnodes(lmn)
                        for k = 1:ndof2
                            rw = ndof2*(connect(a,lmn)-1)+i;
                            cl = ndof2*(connect(b,lmn)-1)+k;
                            K(rw,cl) = K(rw,cl) + k_el(ndof2*(a-1)+i,ndof2*(b-1)+k);
                        end
                    end
                end
            end
        end
    
        % d) Residual vector 
        for a = 1:nelnodes(lmn)
            for i = 1:ndof
                rw = ndof*(connect(a,lmn)-1)+i;
                Res_F(rw,1) = Res_F(rw,1) + Res_F_el(ndof*(a-1)+i,1);
            end
        end

         % e) Internal force vector
        for a = 1:nelnodes(lmn)
            for i = 1:ndof
                rw = ndof*(connect(a,lmn)-1)+i;
                F_int(rw,1) = F_int(rw,1) + f_internal_el(ndof*(a-1)+i,1);
            end
        end

        % f) External force vector
        for a = 1:nelnodes(lmn)
            for i = 1:ndof
                rw = ndof*(connect(a,lmn)-1)+i;
                F_ext(rw,1) = F_ext(rw,1) + f_external_el(ndof*(a-1)+i,1);
            end
        end

       
    
    end


elseif IsProj == 1
% -------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%% PLOTTING PROPERTIES ROUTINE %%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------

    % Create empty entries for solver variables
    J               = [];
    K               = [];
    DRdu            = [];
    Res_F           = [];
    history_var_mat = [];
    F_int           = [];
    F_ext           = [];
    strain_var_mat  = [];

    % Initializing element/nodal properties and gather matrices at zero
    gausspoints_prop_mat = zeros(nelem * func_numberofintegrationpoints,9);
    nodes_prop_mat       = zeros(nnodes,9);
    lmncoord             = zeros(ncoord,maxnodes);
    lmndof               = zeros(ndof,maxnodes);

    % Loop over all the elements
    for lmn = 1:nelem
    
        % Extract coords of nodes, DOF for the current element
        for a = 1:nelnodes(lmn)
            for i = 1:ncoord
                lmncoord(i,a) = coords(i,connect(a,lmn));
            end
            for i = 1:ndof
                lmndof(i,a) = dofs(ndof*(connect(a,lmn)-1)+i);
            end
        end

        % Compute element/nodal properties
        % -----------------------------------------------------------------
        if SolverID == 1
            [~, ~, ~, ~, gausspoints_prop_mat(4*lmn-3:4*lmn,:), nodes_prop_mat_elem, ~, ~, ~] = func_elstif_Local(model_name,lmncoord,lmndof,Delastic,history_var_mat_previousinc(lmn,:),alpha_val,beta_val,e_delta,dmax,n_hood,weights,strain_tolerance,strain_mat_previousinc(lmn,:),1,IsProj,RoutineID);
        % -----------------------------------------------------------------
        elseif SolverID == 2
            [~, ~, ~, gausspoints_prop_mat(4*lmn-3:4*lmn,:), nodes_prop_mat_elem, ~, ~, ~] = func_elstif_Nonlocgradient(model_name,lmncoord,lmndof,Delastic,history_var_mat_previousinc(lmn,:),g,alpha_val,beta_val,e_delta,dmax,strain_tolerance,n_hood,weights,1,IsProj);
        % -----------------------------------------------------------------
        elseif SolverID == 3
            [~, ~, gausspoints_prop_mat(4*lmn-3:4*lmn,:), nodes_prop_mat_elem,~,~] = func_elstif_Nonlocintegral(lmn,lmncoord,lmndof,Delastic,history_var_mat_previousinc(lmn,:),alpha_val,beta_val,e_delta,dmax,n_hood,weights,damage_mat_it,1,IsProj);
        % -----------------------------------------------------------------
        else
            disp("Check your SolverID - globalstiffness")
        end
        
        % Populate the properties matrix
        for a = 1:nelnodes(lmn)
            rw = connect(a,lmn);
            nodes_prop_mat(rw,:) = nodes_prop_mat(rw,:) + nodes_prop_mat_elem(a,:);
        end

    end

    % Calculate the final matrix of nodal properties
    nodes_prop_mat = nodes_prop_mat ./ num_elem_at_node;

else 

    disp("Check your IsProj variable - globalstiffness")

end




end
