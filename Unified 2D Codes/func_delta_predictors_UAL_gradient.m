function [delta_m_bar,delta_u_bar,delta_u_f,delta_f_rct_disp_ebc,delta_e_nl,ArcLength] = func_delta_predictors_UAL_gradient(Jff,Jpp,Jpf,Jfp,Kue_p,Kue_f,Keu_p,Keu_f,Kee,Applied_Displacement_Load,delta_m_bar_0)
% This function uses the predictor values for the deltas as described in
% equation 24, 25

% Define alpha 
alpha=delta_m_bar_0;
in_Jff=inv(Jff);

% Calculate delta predictor values (Equation 24)

% 1. delta_m_bar
delta_m_bar=delta_m_bar_0;

% 2. delta_u_bar
delta_u_bar=alpha*Applied_Displacement_Load;


% 3. delta_u_f
delta_u_f_num=(Kue_f*inv(Kee)*Keu_p*Applied_Displacement_Load*delta_m_bar)-(Jfp*Applied_Displacement_Load*delta_m_bar);
delta_u_f_den=(Jff-(Kue_f*inv(Kee)*Keu_f));

delta_u_f=inv(delta_u_f_den)*delta_u_f_num;

% delta_e_nl
delta_e_nl=(-1)*(inv(Kee))*((Keu_f*delta_u_f)+(Keu_p*Applied_Displacement_Load*delta_m_bar));


%5. delta_f_rct
delta_f_rct_disp_ebc=(-1)*((Jpf*delta_u_f)+(Jpp*Applied_Displacement_Load*delta_m_bar)+(Kue_p*delta_e_nl));

% Calculate Arc Length predictor (Equation 25, spherical only)
ArcLength_0=delta_m_bar*(sqrt((norm(Applied_Displacement_Load,2)))^2+((norm(inv(Jff)*Jpf'*Applied_Displacement_Load))^2));

ArcLength=ArcLength_0;

end