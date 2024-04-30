function [theta_1, theta_2] = func_calc_theta_cosinverse(delta_u_f,del_u_NR,del_m_bar_1,del_m_bar_2,del_u_g,delta_u_bar,Applied_Displacement_Load)
% Calculate the angle using cos inverse 

% Calculate theta_1
m_a1=delta_u_f;
[m_a1_x,m_a1_y] = func_separate_x_y_UAL(m_a1);

m_b1=(delta_u_f+del_u_NR+(del_m_bar_1*del_u_g));
[m_b1_x,m_b1_y] = func_separate_x_y_UAL(m_b1);

m_a2=delta_u_bar;
[m_a2_x,m_a2_y] = func_separate_x_y_UAL(m_a2);

m_b2=(delta_u_bar+(del_m_bar_1*Applied_Displacement_Load));
[m_b2_x,m_b2_y] = func_separate_x_y_UAL(m_b2);

m_num=((m_a1_x'*m_b1_x)+(m_a1_y'*m_b1_y))+((m_a2_x'*m_b2_x)+(m_a2_y'*m_b2_y));
m_den=((norm(m_a1_x,2)*norm(m_b1_x,2))+(norm(m_a1_y,2)*norm(m_b1_y,2)))+((norm(m_a2_x,2)*norm(m_b2_x,2))+(norm(m_a2_y,2)*norm(m_b2_y,2)));      

[t1] = func_cosine_limits_UAL(m_num/m_den);
theta_1=acosd(t1);

% Calculate theta_2
n_a1=delta_u_f;
[n_a1_x,n_a1_y] = func_separate_x_y_UAL(n_a1);

n_b1=(delta_u_f+del_u_NR+(del_m_bar_2*del_u_g));
[n_b1_x,n_b1_y] = func_separate_x_y_UAL(n_b1);

n_a2=delta_u_bar;
[n_a2_x,n_a2_y] = func_separate_x_y_UAL(n_a2);

n_b2=(delta_u_bar+(del_m_bar_2*Applied_Displacement_Load));
[n_b2_x,n_b2_y] = func_separate_x_y_UAL(n_b2);

n_num=((n_a1_x'*n_b1_x)+(n_a1_y'*n_b1_y))+((n_a2_x'*n_b2_x)+(n_a2_y'*n_b2_y));
n_den=((norm(n_a1_x,2)*norm(n_b1_x,2))+(norm(n_a1_y,2)*norm(n_b1_y,2)))+((norm(n_a2_x,2)*norm(n_b2_x,2))+(norm(n_a2_y,2)*norm(n_b2_y,2)));      

[t2] = func_cosine_limits_UAL(n_num/n_den);
theta_2=acosd(t2);

end
