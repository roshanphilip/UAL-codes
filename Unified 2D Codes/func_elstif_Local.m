function [k_el, j_el, damage_el, local_strain_el, gausspoints_prop_mat_elem, nodes_prop_mat_elem, residual_el, f_internal_el, f_external_el] = func_elstif_Local(model_name,lmncoord,dofs,Delastic,damage_mat_previousinc,alpha_val,beta_val,e_delta,dmax,n_hood,weights,strain_tolerance,strain_mat_previousinc,IsM,IsProj,RoutineID)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ============ ELEMENT STIFFNESS MATRIX & INTERNAL FORCE VECTOR ===========
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% -------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SOLVER ROUTINE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
if IsProj == 0

    % Create empty entries for element/nodal properties 
    gausspoints_prop_mat_elem = [];
    nodes_prop_mat_elem = [];   

    % Unroll dofs into a vector
    dofs_vec = dofs(:); % 8x1
    
    % Initialize matrices and vectors at zero
    j_el            = zeros(8,8);
    k_el            = zeros(8,8);
    f_internal_el   = zeros(8,1);
    f_external_el   = zeros(8,1);
    local_strain_el = zeros(1,4);
    damage_el       = zeros(1,4);

    % Setting up: 
    npoints = func_numberofintegrationpoints;   % Number of integration points
    xilist  = func_integrationpoints;           % Positions of integration points
    w       = func_integrationweights;          % Weights of integration points
    
    % Loop over the integration points
    for integ_point = 1:npoints
    
        % ================== SHAPE FUNCTIONS AND DERIVATIVES ==================
        % Compute shape functions && derivatives wrt local coords
        xi = xilist(:,integ_point);             % Size: 2x1
        N = func_shapefunctions(xi);            % Size: 4x1 
        dNdxi = func_shapefunctionderivs(xi);   % Size: 4x2
        
        % Compute the jacobian matrix && its determinant
        dxdxi = lmncoord * dNdxi; % Size: 2x2 - Jacobian matrix
        dt = det(dxdxi);
    
        % Convert shape function derivatives to derivatives wrt global coords
        dNdx = dNdxi/dxdxi; % Size: 4x2 - Multiplies by inverse of dxdxi
    
        % Transpose and expand dNdx (entries of B matrix) to B matrix
        % dNdx is 4x2 ===> B is 3x8
        [B, ~] = func_Bmatrix(dNdx);
    
    
        % ==================== LOCAL STRAIN CALCULATIONS ======================
        % Compute the strain by multiplying the shape function derivatives with
        % the displacements
        strain_vec_gxy = B * dofs_vec; % Size: 3x1 (3x8 * 8x1)
    
        % Extract the strain values for each integration point
        exx = strain_vec_gxy(1,1);
        eyy = strain_vec_gxy(2,1);
        gxy = strain_vec_gxy(3,1);
        
        % Calculate the equivalent strain e_star
        if model_name =='SNS'
            [e_star,s] = func_estar_shear(exx, eyy, gxy);
        else
            [e_star,s] = func_estar(exx, eyy, gxy);
        end
        local_strain_el(1,integ_point) = e_star;

        % Add NaN condition for [s]
        if isnan(s); s(:) = 0; end
    
    
        % =============== MAZAR DAMAGE MODEL IMPLEMENTATION ===============       
        if RoutineID==1     % For UAL only
            % Store the damage variable for each Gauss Point
            if IsM == 0
                e_star_0 = e_star;
                % Calculate damage variable omega
                [omega, domega_dkappa] = func_mazarmodel_Local(e_star_0,alpha_val,beta_val,e_delta,dmax);
                damage_el(1,integ_point)=omega;
            else
                if  ((strain_mat_previousinc(1,integ_point))-(e_star))>strain_tolerance
                    % Apply the Clausius-Duhem inequality
                    e_star_0 = max(strain_mat_previousinc(1,integ_point),e_star);
                    % Calculate damage variable omega
                    [omega, ~] = func_mazarmodel_Local(e_star_0,alpha_val,beta_val,e_delta,dmax);
                    damage_el(1,integ_point)=omega;                   
                    domega_dkappa=0;
                else
                    e_star_0 = e_star;
                    % Calculate damage variable omega
                    [omega, domega_dkappa] = func_mazarmodel_Local(e_star_0,alpha_val,beta_val,e_delta,dmax);
                    damage_el(1,integ_point)=omega;
                end
            end
        elseif RoutineID==2||RoutineID==3    % For NR and FAL            
            % Calculate damage variable omega
            [omega, domega_dkappa] = func_mazarmodel_Local(e_star,alpha_val,beta_val,e_delta,dmax);

            % Store the damage variable for each Gauss Point
            if IsM == 0
                damage_el(1,integ_point) = omega;
            else
                damage_el(1,integ_point) = max(damage_mat_previousinc(1,integ_point),omega);
            end
        end



        % ======================= STRESS CALCULATION ==========================
        stress_vec = (1 - damage_el(1,integ_point)) * Delastic * strain_vec_gxy; % Size: 3x1 (1x1 * 3x3 * 3x1)
    
        
        % ==================== ELEMENT STIFFNESS MATRIX =======================
        % Compute the contribution of this Gauss point to the element stiffness matrix
        KA = (1 - damage_el(1,integ_point)) * B' * Delastic * B;
        
        KB = zeros(8,8);
        for i = 1:4
            for j = 1:2
               row = 2*(i-1)+j;           
               for ii = 1:4
                   for jj = 1:2
                       col = 2 * (ii - 1) + jj;
                       dD_duj = -1 * domega_dkappa * (s(1) * B(1,col) + s(2) * B(2,col) + s(3) * B(3,col));     % dD_duj = - domega_dkappa  * destar/deij * deij/duj
                       for ki = 1:4
                           for kj = 1:2
                               kind = 2 * (ki - 1) + kj;
                               KB(row,col) = KB(row,col) + dD_duj * (1 / (1 - damage_el(1,integ_point))) * KA(row,kind) * dofs_vec(kind);
                           end
                       end
                   end
               end           
            end        
        end
        
        % Assemble element stiffness
        Kint = KA + KB;

        % Compute the contribution of each GP to the element tangent and stiffness matrix 
        j_el = j_el + (Kint* w(integ_point) * dt) ; % Size: 8x8 (1x1 * 8x3 * 3x3 * 3x8 * 1x1 * 1x1)
        k_el = k_el + (KA * w(integ_point) * dt) ; % Size: 8x8 (1x1 * 8x3 * 3x3 * 3x8 * 1x1 * 1x1)
    
        % ================== ELEMENT INTERNAL AND EXTERNAL FORCE VECTOR ====================
        % Compute the contribution of each Gauss Point to the element internal force
        f_internal_el = f_internal_el + B' * stress_vec * w(integ_point) * dt; % Size: 8x1
        f_external_el = zeros(8,1);                                            % Size: 8x1
    
    end

    % ==================== ELEMENT RESIDUAL FORCE VECTOR ======================
    % Assemble the element residual vector.
    residual_el = [f_internal_el - f_external_el];
    


elseif IsProj == 1
% -------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%% PLOTTING PROPERTIES ROUTINE %%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
    
j_el          = [];
k_el          = [];
f_internal_el = [];
f_external_el = [];
residual_el   = [];

    % Unroll dofs into a vector
    dofs_vec = dofs(:); % 8x1
    
    % Initialize matrices and vectors at zero
    local_strain_el           = zeros(1,4);
    damage_el                 = zeros(1,4);
    M                         = zeros(4,4); 
    gausspoints_prop_mat_elem = zeros(4,9);
    prop_gpts                 = zeros(4,9);

    % Setting up: 
    npoints = func_numberofintegrationpoints;   % Number of integration points
    xilist  = func_integrationpoints;           % Positions of integration points
    w       = func_integrationweights;          % Weights of integration points
    
    % Loop over the integration points
    for integ_point = 1:npoints
    
        % ================== SHAPE FUNCTIONS AND DERIVATIVES ==================
        % Compute shape functions && derivatives wrt local coords
        xi = xilist(:,integ_point);             % Size: 2x1
        N = func_shapefunctions(xi);            % Size: 4x1 
        dNdxi = func_shapefunctionderivs(xi);   % Size: 4x2
        
        % Compute the jacobian matrix && its determinant
        dxdxi = lmncoord * dNdxi; % Size: 2x2 - Jacobian matrix
        dt = det(dxdxi);
    
        % Convert shape function derivatives to derivatives wrt global coords
        dNdx = dNdxi/dxdxi; % Size: 4x2 - Multiplies by inverse of dxdxi
    
        % Transpose and expand dNdx (entries of B matrix) to B matrix
        % dNdx is 4x2 ===> B is 3x8
        [B, ~] = func_Bmatrix(dNdx);
    
    
        % ==================== LOCAL STRAIN CALCULATIONS ======================
        % Compute the strain by multiplying the shape function derivatives with
        % the displacements
        strain_vec_gxy = B * dofs_vec; % Size: 3x1 (3x8 * 8x1)
    
        % Extract the strain values for each integration point
        exx = strain_vec_gxy(1,1);
        eyy = strain_vec_gxy(2,1);
        gxy = strain_vec_gxy(3,1);
        exy = gxy/2;

        % Calculate the equivalent strain e_star
        if model_name =='SNS'
            [e_star,~] = func_estar_shear(exx, eyy, gxy);
        else
            [e_star,~] = func_estar(exx, eyy, gxy);
        end
        local_strain_el(1,integ_point) = e_star;


        % =============== MAZAR DAMAGE MODEL IMPLEMENTATION ===============
        % Calculate damage variable omega
        [omega, ~] = func_mazarmodel_Local(e_star,alpha_val,beta_val,e_delta,dmax);
    
        % Store the damage variable for each Gauss Point
        if IsM == 0
            damage_el(1,integ_point) = omega; 
        else
            damage_el(1,integ_point) = max(damage_mat_previousinc(1,integ_point),omega);
        end
        
        % ======================= STRESS CALCULATION ==========================
        stress_vec = (1 - damage_el(1,integ_point)) * Delastic * strain_vec_gxy; % Size: 3x1 (1x1 * 3x3 * 3x1)
    
        % -----------------------------------------------------------------
        % ========== ASSEMBLE PROPERTIES AT ELEMENT GAUSS POINT ===========
        % -----------------------------------------------------------------
        gausspoints_prop_mat_elem(integ_point,:) = [damage_el(1,integ_point) stress_vec' strain_vec_gxy' e_star 0];


        % -----------------------------------------------------------------
        % ===== ASSEMBLE GAUSS POINT PROPERTIES FOR NODAL PROJECTION ======
        % -----------------------------------------------------------------
        for i = 1:4 % number of nodes per element: 4
            prop_gpts(i,1)   = prop_gpts(i,1)   + N(i) * w(integ_point) * dt * damage_el(1,integ_point); % damage
            prop_gpts(i,2:4) = prop_gpts(i,2:4) + N(i) * w(integ_point) * dt * stress_vec';              % stresses [sigma_xx,sigma_yy,tau_xy]
            prop_gpts(i,5:7) = prop_gpts(i,5:7) + N(i) * w(integ_point) * dt * [exx,eyy,exy];            % strains [e_xx, e_yy, gxy]
            prop_gpts(i,8)   = prop_gpts(i,8)   + N(i) * w(integ_point) * dt * e_star;                   % local equivalent strain [e_star]
            prop_gpts(i,9)   = 0;                                                                        % nonlocal equivalent strain             
        end
    
        % Calculate M matrix for projections
        M = M + (N * N' * w(integ_point) * dt); 

    end

    % ====================== COMPUTE NODAL PROPERTIES =====================
    nodes_prop_mat_elem = M \ prop_gpts; 

else

    disp("Check your IsProj variable")

end



end

