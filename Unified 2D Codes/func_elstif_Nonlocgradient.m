function [k_el, j_el, kappa_el, gausspoints_prop_mat_elem, nodes_prop_mat_elem, residual_el, f_internal_el, f_external_el] = func_elstif_Nonlocgradient(model_name,lmncoord,dofs,Delastic,kappa_el_previousinc,g,alpha_val,beta_val,e_delta,dmax,strain_tolerance,n_hood,weights,IsM,IsProj)
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
        
    % Isolate the nodal displacements ux and uy and unroll into a vector
    dofs_u      = dofs(1:2,:);   % Size: 2x4
    dofs_u_vec  = dofs_u(:); % Size: 8x1
    
    % Isolate the nodal nonlocal equivalent strains and unroll into a vector
    dofs_e      = dofs(3,:);     % Size: 1x4
    dofs_e_vec  = dofs_e(:); % Size: 4x1
    
    % Initialize matrices and vectors at zero
    local_strain_el     = zeros(1,4);
    nonlocal_strain_el  = zeros(1,4);
    kappa_el            = zeros(1,4);
    damage_el           = zeros(1,4);
    stress_s1_el        = zeros(1,4);
    f_internal_uu_el    = zeros(8,1);  
    f_external_ee_el    = zeros(4,1);
    k_el                = zeros(8,8);
    j_uu_el             = zeros(8,8);      
    j_eu_el             = zeros(4,8);       
    j_ue_el             = zeros(8,4);        
    j_ee_el             = zeros(4,4);
    
    % Setting up: 
    npoints = func_numberofintegrationpoints;   % Number of integration points
    xilist  = func_integrationpoints;           % Positions of integration points
    w       = func_integrationweights;          % Weights of integration points
    
    % Loop over the integration points
    for integ_point = 1:npoints
    
        % ================== SHAPE FUNCTIONS AND DERIVATIVES ==================
        % Compute shape functions && derivatives wrt local coords. We use the same shape functions for displacement and strain calculations
        xi = xilist(:,integ_point);             % Size: 2x1
        N = func_shapefunctions(xi);            % Size: 4x1 
        dNdxi = func_shapefunctionderivs(xi);   % Size: 4x2
    
        % Compute the jacobian matrix && its determinant
        dxdxi = lmncoord * dNdxi; % Size: 2x2 - Jacobian matrix
        dt = det(dxdxi);
    
        % Convert shape function derivatives to derivatives wrt global coords
        dNdx = dNdxi/dxdxi; % Size: 4x2 
    
        % Construct B matrices from dNdx. Bu (3x8) is used for the displacements 
        % and Be (2x4) is used for the nonlocal equivalent strains
        [Bu, Be] = func_Bmatrix(dNdx);
    
     
        % ==================== LOCAL STRAIN CALCULATIONS ======================
        % Compute the local strain at the Gauss points
        strain_vec_gxy = Bu * dofs_u_vec; % Size: 3x1
        
        % Extract the strain values at each Gauss point
        exx = strain_vec_gxy(1,1);
        eyy = strain_vec_gxy(2,1);
        gxy = strain_vec_gxy(3,1);
        
        % Calculate the equivalent strain e_star and the derivatives {s}
        if model_name =='SNS'
            [e_star,s] = func_estar_shear(exx, eyy, gxy);
        else
            [e_star,s] = func_estar(exx, eyy, gxy);
        end
        local_strain_el(1,integ_point) = e_star;
    
        % Add NaN condition for [s]
        if isnan(s); s(:) = 0; end
        
        
        % ================== NON-LOCAL STRAIN CALCULATIONS ====================
        % Compute and store the nonlocal equivalent strain for each Gauss Point
        nonlocal_strain_el(1,integ_point) = N' * dofs_e_vec; % Size: 1x1
    
    
        % =================== MAZAR MODEL IMPLEMENTATION ======================
        % Calculate the following variables for each Gauss point:
        % a) nonlocal equivalent strain history parameter (kappa)
        % b) damage variable (omega) 
        % c) gradient of omega w.r.t kappa (domega_dkappa)
        [kappa, omega, domega_dkappa] = func_mazarmodel_Nonlocgradient(nonlocal_strain_el(1,integ_point),kappa_el_previousinc(1,integ_point),alpha_val,beta_val,e_delta,dmax,strain_tolerance,IsM);
        
        % Store the output of the Mazar model to the appropriate element vectors
        kappa_el(1,integ_point)  = kappa;
        damage_el(1,integ_point) = omega;
    
        
        % ======================= STRESS CALCULATION ==========================
        stress_vec = (1 - damage_el(1,integ_point)) * Delastic * strain_vec_gxy; % Size: 3x1 (1x1 * 3x3 * 3x1) 
        
        % Calculate the first principal stress and store it to the global stress sigma1 matrix 
        stress_s1_el(1,integ_point) = 0.5*(stress_vec(1,1) + stress_vec(2,1)) + sqrt((0.5*(stress_vec(1,1)-stress_vec(2,1))).^2 + stress_vec(3,1).^2);
        
        
        % ======================== CONSISTENT TANGENT SUBMATRICES =========================
        j_uu_el = j_uu_el + (1 - damage_el(1,integ_point)) * Bu' * Delastic * Bu * w(integ_point) * dt;   % Size: 8x8 (1x1 * 8x3 * 3x3 * 3x8 * 1x1 * 1x1)
        j_ue_el = j_ue_el + Bu' * Delastic * strain_vec_gxy * N' * domega_dkappa * w(integ_point) * dt;   % Size: 8x4 (8x3 * 3x3 * 3x1 * 1x4 * 1x1 * 1x1 * 1x1)
        j_eu_el = j_eu_el + N * s' * Bu * w(integ_point) * dt;                                            % Size: 4x8 (4x1 * 1x3 * 3x8 * 1x1 * 1x1)
        j_ee_el = j_ee_el + (N * N' + g * (Be' * Be)) * w(integ_point) * dt;                              % Size: 4x4 (4x1 * 1x4 + 1x1 * 4x2 * 2x4 * 1x1 * 1x1)
            
        % ======================== STIFFNESS MATRIX =========================
        k_el = j_uu_el;
        
        % ================== ELEMENT INTERNAL FORCE VECTOR ====================
        % Compute the contribution of each Gauss Point to the element internal force 
        f_internal_uu_el = f_internal_uu_el + Bu' * stress_vec * w(integ_point) * dt; % Size: 8x1 (8x3 * 3x1 * 1x1 * 1x1)
        
        % ============ ELEMENT NONLOCAL EXTERNAL FORCE VECTOR =================
        % Compute the contribution of each Gauss Point to the element nonlocal external force vector
        f_external_ee_el = f_external_ee_el + N * e_star * w(integ_point) * dt; % Size: 4x1 (4x1 * 1x1 * 1x1 * 1x1)
    
      
    end
    
    % ========================== STIFFNESS MATRIX =============================
    % Assemble the element stiffness matrix 
    j_el = [j_uu_el -j_ue_el; -j_eu_el j_ee_el];      % Size: 12x12
    
    
    % =============== ELEMENT INTERNAL FORCE VECTOR continued =================
    % Assemble the element internal force 
    f_internal_ee_el = j_ee_el * dofs_e_vec;              % Size: 4x1 (4x4 * 4x1)
    f_internal_el = [f_internal_uu_el; f_internal_ee_el]; % Size: 12x1
    
    
    % =============== ELEMENT EXTERNAL FORCE VECTOR continued =================
    % Assemble the element external force vector. NOTE: The zeros correspond to 
    % f_external_uu, which is the reactions values and they are calculated in 
    % the Newton-Raphson script. With this approach, we contribute with zero 
    % to those values while maintaining the same partition for k_el, 
    % f_internal, f_external. This will allow us to assemble correctly the
    % global external force vector in the globalstiffness script . 
    f_external_el = [zeros(8,1); f_external_ee_el]; 
    
    
    % ==================== ELEMENT RESIDUAL FORCE VECTOR ======================
    % Assemble the element residual vector.
    residual_el = [f_internal_el - f_external_el];
    
    
    % ======== PARTITION STIFFNESS MATRIX AND RESISUAL/FORCES VECTOR ==========
    % Partition stiffness and residual/forces to the element global system:
    % [u1_1, u1_2, u2_1, ... u3_3, u4_3] --> [u1_1, u1_2, u1_3, u2_1, ... u4_3]
    
    % Set up the indices order
    index_series = [1 2 9 3 4 10 5 6 11 7 8 12];
    
    % Partition the matrix and vectors to follow the global system 
    % [u1_1, u1_2, u1_3, u2_1, u2_2, u2_3, ...., unnode_1, unnode_2, unnode_3]
    j_el          = j_el(index_series,index_series);
    residual_el   = residual_el(index_series);




elseif IsProj == 1
% -------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%% PLOTTING PROPERTIES ROUTINE %%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
    
    j_el            = [];
    k_el       = [];
    residual_el     = [];
    kappa_el        = [];
    f_internal_el   = [];
    f_external_el   = [];

    % Isolate the nodal displacements ux and uy and unroll into a vector
    dofs_u      = dofs(1:2,:);   % Size: 2x4
    dofs_u_vec  = dofs_u(:); % Size: 8x1
    
    % Isolate the nodal nonlocal equivalent strains and unroll into a vector
    dofs_e      = dofs(3,:);     % Size: 1x4
    dofs_e_vec  = dofs_e(:); % Size: 4x1
    
    % Initialize matrices and vectors at zero
    local_strain_el           = zeros(1,4);
    nonlocal_strain_el        = zeros(1,4);
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
        % Compute shape functions && derivatives wrt local coords. We use the same shape functions for displacement and strain calculations
        xi = xilist(:,integ_point);             % Size: 2x1
        N = func_shapefunctions(xi);            % Size: 4x1 
        dNdxi = func_shapefunctionderivs(xi);   % Size: 4x2
    
        % Compute the jacobian matrix && its determinant
        dxdxi = lmncoord * dNdxi; % Size: 2x2 - Jacobian matrix
        dt = det(dxdxi);

        % Convert shape function derivatives to derivatives wrt global coords
        dNdx = dNdxi/dxdxi; % Size: 4x2 
    
        % Construct B matrices from dNdx. Bu (3x8) is used for the displacements 
        % and Be (2x4) is used for the nonlocal equivalent strains
        [Bu, ~] = func_Bmatrix(dNdx);
    
     
        % ==================== LOCAL STRAIN CALCULATIONS ======================
        % Compute the local strain at the Gauss points
        strain_vec_gxy = Bu * dofs_u_vec; % Size: 3x1
        
        % Extract the strain values at each Gauss point
        exx = strain_vec_gxy(1,1);
        eyy = strain_vec_gxy(2,1);
        gxy = strain_vec_gxy(3,1);
        exy = gxy/2;
        
        % Calculate the equivalent strain e_star and the derivatives {s}
        if model_name =='SNS'
            [e_star,~] = func_estar_shear(exx, eyy, gxy);
        else
            [e_star,~] = func_estar(exx, eyy, gxy);
        end
        local_strain_el(1,integ_point) = e_star;
    
        
        % ================== NON-LOCAL STRAIN CALCULATIONS ====================
        % Compute and store the nonlocal equivalent strain for each Gauss Point
        nonlocal_strain_el(1,integ_point) = N' * dofs_e_vec; % Size: 1x1
    
    
        % =================== MAZAR MODEL IMPLEMENTATION ======================
        % Calculate the following variables for each Gauss point:
        % a) nonlocal equivalent strain history parameter (kappa)
        % b) damage variable (omega) 
        % c) gradient of omega w.r.t kappa (domega_dkappa)
        [~, omega, ~] = func_mazarmodel_Nonlocgradient(nonlocal_strain_el(1,integ_point),kappa_el_previousinc(1,integ_point),alpha_val,beta_val,e_delta,dmax,strain_tolerance,IsM);
        
        % Store the output of the Mazar model to the appropriate element vectors
        damage_el(1,integ_point) = omega;
    
        
        % ======================= STRESS CALCULATION ==========================
        stress_vec = (1 - damage_el(1,integ_point)) * Delastic * strain_vec_gxy; % Size: 3x1 (1x1 * 3x3 * 3x1) 
        

        % -----------------------------------------------------------------
        % ========== ASSEMBLE PROPERTIES AT ELEMENT GAUSS POINT ===========
        % -----------------------------------------------------------------
        gausspoints_prop_mat_elem(integ_point,:) = [damage_el(1,integ_point) stress_vec' strain_vec_gxy' e_star nonlocal_strain_el(1,integ_point)];


        % -----------------------------------------------------------------
        % ===== ASSEMBLE GAUSS POINT PROPERTIES FOR NODAL PROJECTION ======
        % -----------------------------------------------------------------
        for i = 1:4 % number of nodes per element: 4
            prop_gpts(i,1)   = prop_gpts(i,1)   + N(i) * w(integ_point) * dt * damage_el(1,integ_point);            % damage
            prop_gpts(i,2:4) = prop_gpts(i,2:4) + N(i) * w(integ_point) * dt * stress_vec';                         % stresses [sigma_xx,sigma_yy,tau_xy]
            prop_gpts(i,5:7) = prop_gpts(i,5:7) + N(i) * w(integ_point) * dt * [exx,eyy,exy];                       % strains [e_xx, e_yy, gxy]
            prop_gpts(i,8)   = prop_gpts(i,8)   + N(i) * w(integ_point) * dt * e_star;                              % local equivalent strain [e_star]
            prop_gpts(i,9)   = prop_gpts(i,9)   + N(i) * w(integ_point) * dt * nonlocal_strain_el(1,integ_point);   % nonlocal equivalent strain             
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

