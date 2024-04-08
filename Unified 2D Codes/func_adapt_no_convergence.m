function [flagplot,flaglf,countflaglf,incrflag,loadfactor,dlfactor] = func_adapt_no_convergence(incrflag,loadfactor_stored,dlfactor)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ======== ADAPT LOADFACTOR WHEN NR CONVERGENCE IS NOT CONVERGING =========
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Allow the while loop to continue running in case loadfactor has hit the value 1
flagplot    = 0;
flaglf      = false;
countflaglf = 0;
incrflag    = incrflag + 1;

% Use stored values for load factor (kappa and damage are already being stored)
loadfactor = loadfactor_stored;
        
% Reduce the applied load
dlfactor = dlfactor/2;
loadfactor = min(loadfactor + dlfactor,1);

disp("you just entered highly nonlinear, check results")
disp("========================================")        

end