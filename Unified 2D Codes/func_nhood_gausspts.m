function [neighbours,weights] = func_nhood_gausspts(char_len)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ======= MATRIX WITH NEIGHBORING GAUSS POINTS (NONLOCAL INTEGRAL) ========
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Include global variables
func_include_flags;

% -------------------------------------------------------------------------
lmncoord = zeros(ncoord,maxnodes);
x_coord=zeros(nelem*4,2);


for lmn = 1:nelem
    
    % Extract coords of nodes
    for a = 1:nelnodes(lmn)
        for i = 1:ncoord
            lmncoord(i,a) = coords(i,connect(a,lmn));
        end
    end
    ident = elident_vec(lmn);
    
    % Setting up:
    npoints = func_numberofintegrationpoints;  % Number of integration points
    xilist  = func_integrationpoints;          % Positions of integration points
    
    % Loop over the integration points
    for intpt = 1:npoints
        
        % Compute shape functions && derivatives wrt local coords
        for i = 1:ncoord
            xi(i) = xilist(i,intpt);
        end
        
        N = func_shapefunctions(xi);   % Shape function for each integration point
        x_coord(4 * ident + intpt - 4,:) = N' * lmncoord';
        
    end
    
end

for i = 1:4*nelem
    cnt4nopnear = 0; % counter for no points near
    count = 1;
    for j = 1:4*nelem
        diff_distance = norm(x_coord(i,:) - x_coord(j,:));
        
        if diff_distance <= char_len
            neighbours(i,count) = j;
            weights(i,count)    = exp(-diff_distance / (char_len * char_len));
            count               = count + 1;
            cnt4nopnear         = 1;
        end
        
    end
    
    if cnt4nopnear == 0
        neighbours(i,count) = i;
        weights(i,count)    = 1;
    end
   
    
end

end


