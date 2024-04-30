function [dofs_stored,loadfactor_stored,history_var_mat_stored,incrflag,flagplot,loadfactor,increment,inc_success_counter,flaglf,dlfactor] = func_adapt_slow_convergence(dofs,loadfactor,history_var_mat,last_iteration,dlfactor,increment,inc_success_counter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ======= ADAPT MATRICES AND LOADFACTOR WHEN NR CONVERGENCE IS SLOW =======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Include global variables
func_include_flags;

% Save variables in case of not convergence of next increment
dofs_stored = dofs;
loadfactor_stored = loadfactor;
history_var_mat_stored = history_var_mat;

% Re-initialize at 1 the number of attempts for the increment
incrflag = 1;
flagplot = 1;

% Increment the load
dlfactor            = max_iter/last_iteration * dlfactor;
loadfactor          = min(loadfactor + dlfactor,1);
increment           = increment + 1;
inc_success_counter = inc_success_counter + 1;

if loadfactor == 1
    flaglf = true;
else
    flaglf = false;
end
