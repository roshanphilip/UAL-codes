function [delta_dofs_con, ST, ArcLength_0_FAL, max_ArcLength_FAL, min_ArcLength_FAL, max_displacement_FAL, lambda_con, delta_lambda_con, n_FreeNodes, n_EssentialNodes, reset_counter, u_f_con, delta_f_rct_con, f_rct_con, e_nl_con, delta_u_p_con,max_ArcLength,min_ArcLength,f_reaction_essential_con,m_bar_con,delta_u_bar_con,delta_u_f_con,delta_f_reaction_essential_con,delta_m_bar_con,delta_e_nl_con,ArcLength,n_dofs,c,u_p_con,storage_damage_loadstep_con, Applied_Force_Load, plot_switch, max_delta_lf, min_delta_lf, lf, SolverID,l_Total,n_NodesPerElement,weight_integrationpoint,a,b,n_TotalElements,constraint_type,Applied_Displacement_Load,tolerance,k_max,delta_lf,delta_m_bar_0,DamageThreshholdStrain,n_TotalNodes,n_IntegrationPoints,ConstitutiveMatrix,YoungsModulus,Globalcoordinates_EachNode,connect_LocaltoGlobal_Nodenumbers,n_Global, storage_damage_loadstep, plot_storage, dofs_con, RoutineID, f_reaction] = Inputs_1D()
% This function has the values of all the inputs used in this problem

%=============================================================================================%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% User-defined Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%=============================================================================================%
% RV   - Recommended values
% Model parameters
RoutineID            = 2;                   % 1 - Unified arc length (UAL)   2 - Force-controlled arc length (FAL)  3 - Newton Raphson (NR)
SolverID             = 2;                   % 1 - Local damage   2 - Non-local gradient damage
constraint_type      = 2;                   % 1- Spherical   2- Cylindrical   3- Fixed   
tolerance            = 1e-4;                % Convergence tolerance
ST                   = 1e-8;                % Strain tolerance
lc                   = 6;                   % Characteristic length
k_max                = 50;                  % Max. allowable number of iterations

% Material Parameter
E                    = 30e3;                % Modulus of Elasticity (MPa)
% Geometry parameter
l_Total              = 100.;                % Total Length (mm) of the 1D bar
n_TotalElements      = 151;                 % Total number of elements
damaged_length       = 4;                   % Part of the domain over which the reduction of Young's modulus is applied
notch                = 0.1;                 % Reduction in E of the middle element
notch_type           = 3;                   % 1- No notch    2- Notch in middle element, 3- Notch in middle spread across a fixed length
% Mazar's model parameters
a                    = 0.7;                 % Material parameter for Mazar's damage model
b                    = 1e4;                 % Material parameter for Mazar's damage model
ds                   = 1e-4;                % Damage threshold strain
% Applied load
Applied_Displacement = 0.01;                % (For UAL, NR) Displacement at essential nodes (mm)
Applied_Force_Load   = 5;                   % (For FAL only) Applied load (N)
%=========================================================================%
% UAL parameters
%=========================================================================%
delta_m_bar_0        = 1e-1;                % Initial guess of loadfactor      [RV: 1e-1 to 1e-3]
ArcLength_0          = 1e-4;                % Initial guess of arclength       [RV: 1e-2 to 1e-4]
max_ArcLength        = 1e-4;                % Upper limit of arclength         [RV: < 1e-4]
min_ArcLength        = 1e-24;               % Lower limit of arclength         [RV: 1e-24]
%=========================================================================%
% FAL parameters
% (NOTE: In some cases, especially with local damage, FAL struggles to
% trace the post peak region in a reasonable time frame; Adjust model
% parameters and tolerance accordingly.)
%=========================================================================%
ArcLength_0_FAL      = 1e-4;                % Initial guess of arclength       [RV: 1e-2 to 1e-4]
max_ArcLength_FAL    = 1e-4;                % Upper limit of arclength         [RV: < 1e-4]
min_ArcLength_FAL    = 1e-24;               % Lower limit of arclength         [RV: 1e-24]
max_displacement_FAL = 0.01;                % Max. allowed displacement        (*This is used only as an exit criteria and not in the FAL analysis)
%=========================================================================%
% NR parameters  
% (NOTE: The performance of NR in highly non-linear problems
% is very sensitive to the choice of parameters especially near the peak; If
% convergence is not achieved in the first attempt, rerun the model with
% updated parameters)
%=========================================================================%
delta_lf_0           = 1e-1;                % Initial value of load factor     [RV: 1e-1 to 1e-3]
max_delta_lf         = 1e-4;                % Upper limit of loadfactor        [RV: 1e-4 to 1e-7; choose a smaller value to ensure stability in highly non-linear problems]
min_delta_lf         = 1e-16;               % Lower limit of loadfactor        [RV: 1e-16]
reset_counter        = 5;                   % Number of reset attempts 
%=============================================================================================%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%=============================================================================================%






%=============================================================================================%
%%%%%%%%%%%%%%%%%%%%%%%%% DO NOT CHANGE ANYTHING BEYOND THIS POINT %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%=============================================================================================%






%=============================================================================================%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%=============================================================================================%

Area                 = 1;                                 % Cross-sectional Area (mm2)

n_NodesPerElement    = 2.;                                % Number of nodes in each element

if RoutineID == 1 || RoutineID == 3
    n_EssentialNodes      = 2;                            % Total number of essential nodes UAL and NR routine
elseif RoutineID == 2
    n_EssentialNodes      = 1;                            % Total number of essential nodes for FAL routine
end

weight_integrationpoint   = 2;                            % weight of each integration point

l_EachElement             = l_Total/n_TotalElements;      % Length of each element

% Assign damage threshhold strain based on notch type
if notch_type == 1
    DamageThreshholdStrain = 500;
else
    DamageThreshholdStrain = ds;
end

delta_lf                  = delta_lf_0;
if SolverID == 1    % For UAL
    ArcLength                 = ArcLength_0;
elseif SolverID == 2
    ArcLength                 = ArcLength_0_FAL;
end

% Calculate non-local gradient parameter
c                         = (lc^2)/2;

Applied_Displacement_Load = [0; Applied_Displacement];
n_TotalNodes              = ((n_NodesPerElement-1) * n_TotalElements) + 1; % Total number of nodes
n_IntegrationPoints       = n_NodesPerElement - 1;                         % Number of integration points per element
n_FreeNodes               = n_TotalNodes - n_EssentialNodes;               % Total number of free nodes

% Calculate number of dofs and initialize dofs vector
if SolverID == 1 % Local damage
    n_dofs                = 1;                                             % No. of degrees of freedom per node
    [dofs_con, u]         = deal(zeros((n_dofs * n_TotalNodes),1));
elseif SolverID == 2 % Gradient damage
    n_dofs                = 2;                                             % No. of degrees of freedom per node
    [dofs_con]            = deal(zeros(((n_dofs) * n_TotalNodes),1));
    [u]                   = deal(zeros(((n_dofs-1) * n_TotalNodes),1));
end

u_p_con               = zeros(n_EssentialNodes,1);
u_f_con               = zeros(n_FreeNodes,1);
delta_dofs_con        = zeros((n_dofs * n_TotalNodes),1);
[ConstitutiveMatrix,YoungsModulus] = func_ConstitutiveMatrix(n_TotalElements,E,Area,l_EachElement,damaged_length,l_Total,notch_type,notch);
[lambda_con, delta_lambda_con, delta_f_rct_con, f_rct_con, e_nl_con, delta_u_p_con, f_reaction_essential_con,m_bar_con,delta_u_bar_con,delta_u_f_con,delta_f_reaction_essential_con,delta_m_bar_con,delta_e_nl_con] = deal(zeros());

% Defining the Global coordinates of each node in the +ve X direction
Globalcoordinates_EachNode = zeros(1,n_TotalNodes);
for i = 2 : 1: n_TotalNodes
    Globalcoordinates_EachNode(1,i) = Globalcoordinates_EachNode(1,i-1) + l_EachElement;
end

% Connecting local and global node numbers
connect_LocaltoGlobal_Nodenumbers = zeros(n_NodesPerElement,n_TotalElements);
for i = 1 : 1 : n_TotalElements
    if  n_NodesPerElement == 2
        connect_LocaltoGlobal_Nodenumbers(1,i) = i;
        connect_LocaltoGlobal_Nodenumbers(2,i) = i + 1;
    end
end

% Global unit normal vector matrix
n_Global                                = zeros(n_EssentialNodes,n_TotalNodes);
n_Global(1,1)                           = 1;
n_Global(n_EssentialNodes,n_TotalNodes) = -1;

% Initializing variables
[f_reaction,lf]                         = deal(0);

% Storage matrix for the damage, strain and non local strain in each element across the 1D bar within a loadstep, in each element across the 1D bar
% at all converged loadsteps and all converged loadstep values respectively
[storage_damage_loadstep_con,storage_damage_loadstep]=deal(zeros(1,n_TotalElements));

% Initialize plot variables
plot_storage                            = zeros(1,2);
plot_switch                             = 1;

end
