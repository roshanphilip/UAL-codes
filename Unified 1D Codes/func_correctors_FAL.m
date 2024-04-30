function [del_lamda,imaginary_roots] = func_correctors_FAL(Jff,delta_lamda_lastconverged,delta_u_p_con,delta_u_f_con,delta_e_nl_con,pos_essential_dofs,del_dofs_A_free,del_dofs_B_free, q_f,Beta,delta_lamda,delta_u,Arc_Length,step5_counter,ID_prescribed_u,ID_free)

%Unroll delta_u
delta_u_EssentialNodes = delta_u(ID_prescribed_u,1);
delta_u_FreeNodes      = delta_u(ID_free,1);

Alpha1=(del_dofs_B_free'*del_dofs_B_free)+((Beta^2)*q_f'*q_f);
Alpha2=(2*(Beta^2)*delta_lamda*q_f'*q_f)+(del_dofs_A_free'*del_dofs_B_free)+(del_dofs_B_free'*del_dofs_A_free)+(delta_u_FreeNodes'*del_dofs_B_free)+(del_dofs_B_free'*delta_u_FreeNodes);
Alpha3=(delta_u_FreeNodes'*delta_u_FreeNodes)+((Beta^2)*(delta_lamda^2)*q_f'*q_f)-(Arc_Length^2)+(del_dofs_A_free'*del_dofs_A_free)+(delta_u_FreeNodes'*del_dofs_A_free)+(del_dofs_A_free'*delta_u_FreeNodes);

% Solve the quadratic equation
del_lamda_1=((-1*Alpha2)+sqrt((Alpha2^2)-(4*Alpha1*Alpha3)))/(2*Alpha1);
del_lamda_2=((-1*Alpha2)-sqrt((Alpha2^2)-(4*Alpha1*Alpha3)))/(2*Alpha1);

% Calculate the two del_u values
del_u_free_1=del_dofs_A_free+(del_dofs_B_free*del_lamda_1);
del_u_free_2=del_dofs_A_free+(del_dofs_B_free*del_lamda_2);

% Calculate the two P values
P1=((delta_u_FreeNodes+del_u_free_1)'*delta_u_f_con)+((Beta^2)*delta_lamda_lastconverged*(delta_lamda+del_lamda_1)*q_f'*q_f);
P2=((delta_u_FreeNodes+del_u_free_2)'*delta_u_f_con)+((Beta^2)*delta_lamda_lastconverged*(delta_lamda+del_lamda_2)*q_f'*q_f);

% Calculate the two DOT values
DOT1=((delta_u_FreeNodes+del_u_free_1)'*delta_u_FreeNodes)+((Beta^2)*delta_lamda*(delta_lamda+del_lamda_1)*q_f'*q_f);
DOT2=((delta_u_FreeNodes+del_u_free_2)'*delta_u_FreeNodes)+((Beta^2)*delta_lamda*(delta_lamda+del_lamda_2)*q_f'*q_f);

if step5_counter==1
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



if isreal(del_lamda)~=1 || isreal(P1)~=1 ||isreal(P2)~=1
    imaginary_roots=1;
else
    imaginary_roots=0;
end

end

