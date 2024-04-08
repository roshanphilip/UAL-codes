function [del_lamda,imaginary_roots] = func_solve_arclength_eq_FAL(Jff,delta_lamda_lastconverged,delta_u_last_converged,del_u_A_free,del_u_B_free,q_f,Beta,delta_lamda,delta_u,ArcLength,counter,ID_dofs_list_at_ebc,ID_dofs_list_at_nbc,increment,ArcLength_0)
% This function solves the arc length equation 

% Set ArcLength to its initial guess value at the first increment 
if increment == 1 
    ArcLength = ArcLength_0;
end

% Unroll delta_u
delta_u_FreeNodes = delta_u(ID_dofs_list_at_nbc,1);

% Unroll delta_u_FreeNodes
delta_u_last_converged_free = delta_u_last_converged(ID_dofs_list_at_nbc,1);

% Calculate alpha1, alpha2 and alpha3 coefficients of the arc length equation 
Alpha1=(del_u_B_free'*del_u_B_free)+((Beta^2)*q_f'*q_f);
Alpha2=(2*(Beta^2)*delta_lamda*q_f'*q_f)+(del_u_A_free'*del_u_B_free)+(del_u_B_free'*del_u_A_free)+(delta_u_FreeNodes'*del_u_B_free)+(del_u_B_free'*delta_u_FreeNodes);
Alpha3=(delta_u_FreeNodes'*delta_u_FreeNodes)+((Beta^2)*(delta_lamda^2)*q_f'*q_f)-(ArcLength^2)+(del_u_A_free'*del_u_A_free)+(delta_u_FreeNodes'*del_u_A_free)+(del_u_A_free'*delta_u_FreeNodes);

% Calcuate the roots of the quadratic equation
del_lamda_1=((-1*Alpha2)+sqrt((Alpha2^2)-(4*Alpha1*Alpha3)))/(2*Alpha1);
del_lamda_2=((-1*Alpha2)-sqrt((Alpha2^2)-(4*Alpha1*Alpha3)))/(2*Alpha1);

%=========================================================================%
%============= Parameters used in choosing the correct roots =============%
% Ref: Michael A Crisfield. A fast incremental/iterative solution procedure that handles "snap-through”. In Computational methods in nonlinear structural and solid mechanics, pages 55–62. Elsevier, 1981.
%=========================================================================%
% Calculate the two del_u values
del_u_free_1=del_u_A_free+(del_u_B_free*del_lamda_1);
del_u_free_2=del_u_A_free+(del_u_B_free*del_lamda_2);

% Calculate the two P values
P1=((delta_u_FreeNodes+del_u_free_1)'*delta_u_last_converged_free)+((Beta^2)*delta_lamda_lastconverged*(delta_lamda+del_lamda_1)*q_f'*q_f);
P2=((delta_u_FreeNodes+del_u_free_2)'*delta_u_last_converged_free)+((Beta^2)*delta_lamda_lastconverged*(delta_lamda+del_lamda_2)*q_f'*q_f);

% Calculate the two DOT values
DOT1=((delta_u_FreeNodes+del_u_free_1)'*delta_u_FreeNodes)+((Beta^2)*delta_lamda*(delta_lamda+del_lamda_1)*q_f'*q_f);
DOT2=((delta_u_FreeNodes+del_u_free_2)'*delta_u_FreeNodes)+((Beta^2)*delta_lamda*(delta_lamda+del_lamda_2)*q_f'*q_f);
%=========================================================================%

% Choose the correct root of the quadratic equation 
if counter==1
    sign_det=det(Jff);

    if sign(sign_det)==sign(del_lamda_1)
        del_lamda=del_lamda_1;
    else
        del_lamda=del_lamda_2;
    end

else
    if DOT1>DOT2
        del_lamda=del_lamda_1;
    elseif DOT2>DOT1
        del_lamda=del_lamda_2;
    elseif DOT1==DOT2
        if del_lamda_1>0 && del_lamda_2<0
            del_lamda=del_lamda_1;
        else
            del_lamda=del_lamda_2;
        end
    else
        del_lamda=del_lamda_1;
    end
end

% Check for imaginary roots 
if isreal(del_lamda)~=1 || isreal(P1)~=1 ||isreal(P2)~=1
    imaginary_roots=1;
else
    imaginary_roots=0;
end

end

