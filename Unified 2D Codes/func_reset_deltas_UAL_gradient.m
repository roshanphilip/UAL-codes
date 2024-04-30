function [delta_m_bar,delta_u_bar,delta_u_f,delta_f_rct_disp_ebc,delta_e_nl] = func_reset_deltas_UAL_gradient(convergance_flag,delta_m_bar,delta_m_bar_conv,delta_u_bar_conv,delta_u_f_conv,delta_f_rct_disp_ebc_conv,delta_e_nl_conv,Applied_Displacement_Load)
% This function resets the values of all the deltas based on the current state of the system  
    
% A. If converged in the previous loadstep, set the deltas to the last converged values of deltas from the previous loadstep  
if convergance_flag==1 
    delta_m_bar              = delta_m_bar_conv;
    delta_u_f                = delta_u_f_conv;
    delta_u_bar              = delta_u_bar_conv;
    delta_f_rct_disp_ebc     = delta_f_rct_disp_ebc_conv;
    delta_e_nl               = delta_e_nl_conv;
    
% B. If convergence not reached in the last loadstep after k_max iterations, reduce the last converged deltas by a factor of x (obtained from eq 37)
elseif convergance_flag==0 
    delta_m_bar              = delta_m_bar_conv;
    delta_u_bar              = delta_m_bar*Applied_Displacement_Load;
    delta_u_f                = delta_u_f_conv;
    delta_f_rct_disp_ebc     = delta_f_rct_disp_ebc_conv;
    delta_e_nl               = delta_e_nl_conv;
end

end