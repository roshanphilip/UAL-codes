function [delta_m_bar,delta_u_bar,delta_u_f,delta_f_reaction_essential,delta_e_nl] = func_delta_predictors_UAL(n, Applied_Displacement_Load,delta_m_bar_0,delta_m_bar_con,J,SolverID,ID_prescribed_u,ID_free_u,ID_nl_strain)
% This function calculates the predictor value for UAL analysis

if SolverID == 1 % Local Damage
    delta_e_nl = []; 

    % Partition J
    Jee                            =   J(ID_prescribed_u,ID_prescribed_u);
    Jef                            =   J(ID_prescribed_u,ID_free_u);
    Jfe                            =   J(ID_free_u,ID_prescribed_u);
    Jff                            =   J(ID_free_u,ID_free_u);
    
    % Calculate delta predictor values 
    if n == 1 
        alpha                      =   delta_m_bar_0;
    else                        
        alpha                      =   delta_m_bar_con;
    end

    delta_m_bar                    =   alpha;
    delta_u_bar                    =   delta_m_bar*Applied_Displacement_Load;
    delta_u_f                      =   -1*delta_m_bar_0*(Jff\Jfe)*Applied_Displacement_Load;
    delta_f_reaction_essential     =   delta_m_bar_0*((Jef*(Jff\Jfe))-Jee)*Applied_Displacement_Load;
    

elseif SolverID == 2 % Non-local Gradient 

    % Partition J
    Jee                            =   J(ID_nl_strain,ID_nl_strain);
    Jeu_p                          =   J(ID_nl_strain,ID_prescribed_u);
    Jue_p                          =   J(ID_prescribed_u,ID_nl_strain);
    Jue_f                          =   J(ID_free_u,ID_nl_strain);
    Jeu_f                          =   J(ID_nl_strain,ID_free_u);
                   
    JUff                           =   J(ID_free_u,ID_free_u);
    JUfp                           =   J(ID_free_u,ID_prescribed_u);
    JUpf                           =   J(ID_prescribed_u,ID_free_u);
    JUpp                           =   J(ID_prescribed_u,ID_prescribed_u);
   
    % Calculate delta predictor values 
    if n == 1 
        alpha                      =   delta_m_bar_0;
    else                        
        alpha                      =   delta_m_bar_con;
    end

    delta_m_bar                    =   alpha;
    delta_u_bar                    =   delta_m_bar_0 * Applied_Displacement_Load;
                   
    delta_u_f_num                  =   (Jue_f * (Jee\Jeu_p) * Applied_Displacement_Load * delta_m_bar) - (JUfp * Applied_Displacement_Load * delta_m_bar);
    delta_u_f_den                  =   (JUff -(Jue_f * (Jee\Jeu_f)));
    delta_u_f                      =   delta_u_f_den\delta_u_f_num;
                   
    delta_e_nl                     =   -Jee\((Jeu_f*delta_u_f)+(Jeu_p*Applied_Displacement_Load*delta_m_bar));
    
    delta_f_reaction_essential     =   (-1)*((JUpf*delta_u_f)+(JUpp*Applied_Displacement_Load*delta_m_bar)+(Jue_p*delta_e_nl));

end
end