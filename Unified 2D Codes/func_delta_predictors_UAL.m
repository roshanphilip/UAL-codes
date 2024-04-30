function [delta_m_bar,delta_u_bar,delta_u_f,delta_f_reaction_essential,ArcLength_0,ArcLength] = func_delta_predictors_UAL(Jee,Jef,Jfe,Jff,Applied_Displacement_Load,Beta,del_u_g,del_u_NR,n,r_e,r_f,g,delta_u_bar_con,delta_u_f_con,delta_f_reaction_essential_con,delta_m_bar_con,ArcLength_0,ArcLength)
% This function uses the predictor values for the deltas as described in
% equation 24, 25

if n==1
        
    ArcLength=ArcLength_0;
    
    % Calculate del_m_bar
    x1=(norm(Applied_Displacement_Load,2))^2;
    x2=(norm(del_u_g,2))^2;
    x3=(Jee*Applied_Displacement_Load)+(Jef*del_u_g);
    x4=(Beta^2)*((norm(x3,2))^2);

    del_m_bar_1=ArcLength/(sqrt(x1+x2+x4));
    del_m_bar_2=-(ArcLength/(sqrt(x1+x2+x4)));

    del_m_bar=max(del_m_bar_1,del_m_bar_2);

else   

    % Calculate del_m_bar
    x1=(norm(Applied_Displacement_Load,2))^2;
    x2=(norm(del_u_g,2))^2;
    x3=(Jee*Applied_Displacement_Load)+(Jef*del_u_g);
    x4=(Beta^2)*((norm(x3,2))^2);

    del_m_bar_1=ArcLength/(sqrt(x1+x2+x4));
    del_m_bar_2=-(ArcLength/(sqrt(x1+x2+x4)));
    
    % Identify sign of delta_m_bar
    sign=(delta_u_f_con'*del_u_g)+(delta_u_bar_con'*Applied_Displacement_Load);

    if sign>0
        del_m_bar=max(del_m_bar_1,del_m_bar_2);
    elseif sign<0
        del_m_bar=min(del_m_bar_1,del_m_bar_2);
    end

end

% del_u_bar=del_m_bar*Applied_Displacement_Load;
del_u_f=del_u_NR+(del_m_bar*del_u_g);
del_f_reaction=-(r_e+(del_m_bar*Jee*Applied_Displacement_Load)+(Jef*del_u_f));

delta_m_bar=del_m_bar;
delta_u_bar=delta_m_bar*Applied_Displacement_Load;
delta_u_f=del_u_f;
delta_f_reaction_essential=del_f_reaction;

end