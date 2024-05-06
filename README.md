# UAL-codes
The repository contains the MATLAB code used to model the examples in "A new unified arc-length method for damage mechanics problems" [doi:10.1007/s00466-024-02473-5]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ref: R.P. Saji, P. Pantidis, M. Mobasher, A new unified arc-length method for damage mechanics problems 
% Journal: Computational Mechanics [doi:10.1007/s00466-024-02473-5]
% Author 1: Roshan Philip Saji
% email: rs7625@nyu.edu
% Author 2: Dr. Panos Pantidis
% email: pp2624@nyu.edu
% Author 3: Prof. Mostafa Mobasher 
% email: mostafa.mobasher@nyu.edu
%
% Computational Solid Mechanics Lab, New York University Abu Dhabi
% Date: 26-Feb-2024 
% Link to manuscript: https://doi.org/10.48550/arXiv.2308.13758
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This is the Read me/License file for the 1D and 2D numerical examples presented in the above referenced manuscript.
 
Here, we present a novel unified arc-length method (UAL) for continuum damage mechanics problems. UAL is able to combine the advantages of the commonly used 
Newton-Raphson (NR) solver and the traditional force-controlled arc-length (FAL) solver that make a trade-off between computational cost and solution accuracy when 
solving problems with complex equilibrium paths. The UAL solver has proven to be more robust than the NR and FAL solvers in terms of the computational time taken and its ability
to capture complex snap-back and snap-through regions in the equilibrium path.

This code package contains the following:
• 1D : Solvers     - UAL, FAL, NR 
       Damage laws - Local, Non-Local Gradient
• 2D : Solvers     - UAL, NR
       Damage laws - Local, Non-Local Gradient
• A figure generator for creating crossplots for various scenarios  

We invite you to try out different continuum damage mechanics problems using this package and experience first-hand the benefits of our novel UAL solver.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%================================================================================================================================================================================%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%================================================================================================================================================================================%
I] Steps to run 1D codes:

1. Open the folder named "Unified 1D Codes"
2. Open input file - "Inputs_1D.m"
    • Update the user-defined parameters based on the values in Tables 1 and 3 of the manuscript [Valid for examples presented in Section 6.1 of manuscript]
3. Open main script file - "FEM_Main_Script_1D.m"
	• Enter the file path and results storage folder path
	• Enter results file name 
	• Run main script file - "FEM_Main_Script_1D.m"
	
%================================================================================================================================================================================%
II] Steps to run 2D codes:

1. Open the folder named "Unified 2D Codes"
2. Open main script file - "FEM_Main_Script_2D.m"
	• Enter the file path and results storage folder path
	• Enter the model name 
	  ♦ SNT - Single Notch Tension 
      ♦ SSNT - Symmetric Single Notch Tension (Coarse and Fine)
      ♦ TNT - Two Notch Tension (Coarse and Fine)
      ♦ SNS - Single Notch Shear
    • Enter the following model parameters mentioned in Section 6.2. (refer subsection 6.2.1 to 6.2.4 for model specific parameters)
	  - Convergence tolerance
	  - ArcLength_upper_limit 
	  - ArcLength_lower_limit 
	  - Scheme ID 
	  - Strain tolerance 
3. Run main script file - "FEM_Main_Script_2D.m"

%================================================================================================================================================================================%
III] Steps to Plotting function code:

1. Open the folder named "Plotting functions"
2. Open the file - "Figure_gererator.m"
   • Enter the file path
   • Enter the figure name 
   • Enter the axis labels
   • Enter the axis limits
   • Choose one of the following plot types and comment out the rest:
      A. Plot with markers 
      B. Plot without markers 
      C. Self comparison plots (*useful when comparing results from the same model under different input parameters)
   • Choose the features needed from the following and comment out the rest:
      X1. Zoom into a portion of the graph
	  X2. Text marker annotation
	  X3. Legends 
	  X4. Scientific notation 
3. Run the file - "Figure_gererator.m"
*Ensure that the data files are stored in the same folder as "Figure_gererator.m"
%================================================================================================================================================================================%
NOTE:
* Ensure that you are using the latest version of MATLAB (R2023b or later versions).
* Rename the variable 'connect' to 'connect_xyz' throughout the code if you don't have access to the latest version of MATLAB.
* Ensure that you've created a folder to save the data files at each increment.  
%================================================================================================================================================================================%


