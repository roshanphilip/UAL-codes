function [damage_mat] = func_localdamagemat(dofs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ==================== NEWTON RAPHSON ANALYSIS METHOD =====================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Include global variables
func_include_flags;

% -------------------------------------------------------------------------
damage_mat = zeros(nelem,4);
lmncoord   = zeros(ncoord,maxnodes);
lmndof     = zeros(ndof,maxnodes);

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
    
    % Setting up:
    npoints    = func_numberofintegrationpoints;  % Number of integration points
    xilist     = func_integrationpoints;          % Positions of integration points
    lmndof_vec = lmndof(:);
    ident      = elident_vec(lmn);
    
    % Loop over the integration points
    for intpt = 1:npoints
        
        % Compute shape functions && derivatives wrt local coords
        xi = xilist(:,intpt); 
        dNdxi = func_shapefunctionderivs(xi); % Derivative of shape function wrt to local coordinate for each integration point
        
        % Compute the jacobian matrix && its determinant
        dxdxi = lmncoord * dNdxi;
       
        % Convert shape function derivatives to derivatives wrt global coords
        dNdx = dNdxi/dxdxi;
        
        % Transpose and expand dNdx (entries of B matrix) to B matrix 
        % dNdx is 4x2 ===> B is 3x8
        [B, ~] = func_Bmatrix(dNdx);
    
        % ==================== LOCAL STRAIN CALCULATIONS ======================
        % Compute the infinitesimal strain by multiplying the displacements
        % with the shape function derivatives matrix
        strain_vec_gxy = B * lmndof_vec; % Size: 3x1 (3x8 * 8x1)       
        
        % Extract the strain values for each integration point
        exx = strain_vec_gxy(1,1);
        eyy = strain_vec_gxy(2,1);
        gxy = strain_vec_gxy(3,1);
        
        % Calculate the equivalent strain e_star
        [e_star,~] = func_estar(exx, eyy, gxy);


        % =============== MAZAR DAMAGE MODEL IMPLEMENTATION ===============
        % Calculate damage variable omega
        [omega, ~] = func_mazarmodel_Nonlocintegral(e_star,alpha_val,beta_val,e_delta,dmax);

        % Store the damage variable for each Gauss Point
        damage_mat(ident,intpt) = max(damage_mat(ident,intpt),omega);    
        
    end

end