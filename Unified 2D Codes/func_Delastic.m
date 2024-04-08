function [Delastic] = func_Delastic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =============== CONSTITUTIVE MATRIX (3x3 for 2D problems) ===============
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Include global variables
func_include_flags;

% Extract material properties
mu = materialprops(1);          % Shear modulus
nu = materialprops(2);          % Poisson's ratio
E = 2*mu*(1+nu);                % Young's modulus
planestrain = materialprops(3); % Planestress vs planestrain

% Note: planestrain = 0 => plane stress 
%       planestrain = 1 => plane strain

% Initialize the constitutive matrix D with zeros (3x3 for 2D problems)
Delastic = zeros(3,3);

% Compute the non-zero components of the elastic constitutive matrix D 
if planestrain == 0     % plane stress formulation

    Delastic(1,1) = E/(1-nu^2);
    Delastic(2,2) = E/(1-nu^2);
    Delastic(3,3) = E/(1-nu^2) * (1-nu)/2;
    Delastic(1,2) = E/(1-nu^2) * nu;
    Delastic(2,1) = E/(1-nu^2) * nu;

elseif planestrain == 1 % plane strain formulation

    Delastic(1,1) = E/((1+nu)*(1-2*nu)) * (1-nu);
    Delastic(2,2) = E/((1+nu)*(1-2*nu)) * (1-nu);
    Delastic(3,3) = E/((1+nu)*(1-2*nu)) * (1-2*nu)/2;
    Delastic(1,2) = E/((1+nu)*(1-2*nu)) * nu;
    Delastic(2,1) = E/((1+nu)*(1-2*nu)) * nu;

else % neither plane stress nor plane strain 

    print("Check your inputs in the material properties")

end