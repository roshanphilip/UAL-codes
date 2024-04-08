function [del_f_reaction, del_u_f, del_m_bar,exit_counter] = func_solve_system_equations_nonconsistent_partitioned(delta_u_bar,delta_u_f,delta_f_reaction_essential,J_F,J_E,J_EF,J_FE,Res_F_E_rct,Res_F_F,Beta,Applied_Displacement_Load,ArcLength)
% This function finds the values of corrections for f_reaction, u_f and m_bar

% Initialize imaginary roots counter
exit_counter=0;

% Calculate del_m_bar
del_u_NR=(-1)*J_F\Res_F_F;
del_u_g=(-1)*(J_F)\(J_FE*Applied_Displacement_Load);

% Calculate a
a1=(J_E*Applied_Displacement_Load)+(J_EF*del_u_g);
a=((norm(Applied_Displacement_Load,2))^2)+((norm(del_u_g,2))^2)+((Beta^2)*((norm(a1,2))^2));

% Calculate b
b1=delta_u_bar'*Applied_Displacement_Load;
b2=(delta_u_f+del_u_NR)'*del_u_g;
b3=b1+b2;
b4=(delta_f_reaction_essential-Res_F_E_rct-(J_EF*del_u_NR));
b5=(J_E*Applied_Displacement_Load)+(J_EF*del_u_g);
b6=(Beta^2)*b4'*b5;
b=2*(b3+b6);

% Calculate c
c1=(norm(delta_u_bar,2))^2;
c2=(norm((delta_u_f+del_u_NR),2))^2;
c3=(Beta^2)*((norm((delta_f_reaction_essential-Res_F_E_rct-(J_EF*del_u_NR)),2))^2);
c=c1+c2+c3-(ArcLength^2);

% Calculate d
d=(b^2-4*(a*c));

%========== Check for imaginary roots ==========%
if d<0
    del_f_reaction=[];
    del_u_f=[];
    del_m_bar=[];
    exit_counter=1;
    disp("IMAGINARY ROOTS DETECTED")
    disp("========================================")
    disp("CONSIDER CHANGING CONTROL PARAMETERS LIKE ARCLENGTH LIMIT, STRAIN TOLERANCE OR MAXIMUM ALLOWABLE DAMAGE")
    disp("========================================")
    return
end

% Calculate the two del_m_bar values
del_m_bar_1=(-b+sqrt(d))/(2*a);
del_m_bar_2=(-b-sqrt(d))/(2*a);

% Calculate the two angles using cos inverse
[theta_1, theta_2] = func_calc_theta_cosinverse(delta_u_f,del_u_NR,del_m_bar_1,del_m_bar_2,del_u_g,delta_u_bar,Applied_Displacement_Load);

% Choose the smaller angle
if theta_1<theta_2
    del_m_bar=del_m_bar_1;
else
    del_m_bar=del_m_bar_2;
end

% Calculate corrector value of u_f
del_u_f=del_u_NR+(del_m_bar*del_u_g);

% Calculate corrector value of f_reaction 
del_f_reaction=(-1)*(Res_F_E_rct+(del_m_bar*J_E*Applied_Displacement_Load)+(J_EF*del_u_f));

end