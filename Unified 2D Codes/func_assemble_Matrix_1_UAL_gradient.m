function [Matrix_1] = func_assemble_Matrix_1_UAL_gradient(M_bar_u,M_u,M_bar_e,M_e,sub_1,sub_2,sub_3,sub_4,sub_5,sub_6,sub_7,sub_8,sub_9,sub_10,sub_11,sub_12,sub_13,sub_14,sub_15,sub_16)
% This function assembles the Matrix_1 

% % Initialize flags matrix
flags=zeros((M_bar_u+M_u+M_bar_e+M_e+1),1);

% % Initialize Matrix 1
Matrix_1=zeros((M_bar_u+M_u+M_bar_e+M_e+1),(M_bar_u+M_u+M_bar_e+M_e+1));

% % Assign an identifier to the position of component matrices within the intended assembled matrix
for i=1:1:M_bar_u
    % For sub matrix 1
    for j=1:1:M_bar_u
        flags(i,j)=1;
    end
    % For sub matrix 2
    for j=M_bar_u+1:1:M_bar_u+M_u
        flags(i,j)=2;
    end
    % For sub matrix 3
    for j=M_bar_u+M_u+1        
        flags(i,j)=3;
    end
    % For sub matrix 4
    for j=M_bar_u+M_u+2:1:(M_bar_u+M_u+M_bar_e+M_e+1)     
        flags(i,j)=4;
    end
end


for i=M_bar_u+1:1:M_bar_u+M_u
    % For sub matrix 5
    for j=1:1:M_bar_u
        flags(i,j)=5;
    end
    % For sub matrix 6
    for j=M_bar_u+1:1:M_bar_u+M_u
        flags(i,j)=6;
    end
    % For sub matrix 7
    for j=M_bar_u+M_u+1        
        flags(i,j)=7;
    end
    % For sub matrix 8
    for j=M_bar_u+M_u+2:1:(M_bar_u+M_u+M_bar_e+M_e+1)       
        flags(i,j)=8;
    end
end


for i=M_bar_u+M_u+1
    % For sub matrix 9
    for j=1:1:M_bar_u
        flags(i,j)=9;
    end
    % For sub matrix 10
    for j=M_bar_u+1:1:M_bar_u+M_u
        flags(i,j)=10;
    end
    % For sub matrix 11
    for j=M_bar_u+M_u+1        
        flags(i,j)=11;
    end
    % For sub matrix 12
    for j=M_bar_u+M_u+2:1:(M_bar_u+M_u+M_bar_e+M_e+1)     
        flags(i,j)=12;
    end
end

for i=M_bar_u+M_u+2:1:(M_bar_u+M_u+M_bar_e+M_e+1)  
    % For sub matrix 13
    for j=1:1:M_bar_u
        flags(i,j)=13;
    end
    % For sub matrix 14
    for j=M_bar_u+1:1:M_bar_u+M_u
        flags(i,j)=14;
    end
    % For sub matrix 15
    for j=M_bar_u+M_u+1        
        flags(i,j)=15;
    end
    % For sub matrix 16
    for j=M_bar_u+M_u+2:1:(M_bar_u+M_u+M_bar_e+M_e+1)     
        flags(i,j)=16;
    end
end


% Identify the position of submatrices 1 to 9 within the flags matrix
[Ar,Ac]=find(flags==1);
[Br,Bc]=find(flags==2);
[Cr,Cc]=find(flags==3);
[Dr,Dc]=find(flags==4);
[Er,Ec]=find(flags==5);
[Fr,Fc]=find(flags==6);
[Gr,Gc]=find(flags==7);
[Hr,Hc]=find(flags==8);
[Ir,Ic]=find(flags==9);
[Jr,Jc]=find(flags==10);
[Kr,Kc]=find(flags==11);
[Lr,Lc]=find(flags==12);
[Mr,Mc]=find(flags==13);
[Nr,Nc]=find(flags==14);
[Or,Oc]=find(flags==15);
[Pr,Pc]=find(flags==16);


% Assemble Matrix_1
Matrix_1(Ar(1):Ar(size(Ar)),Ac(1):Ac(size(Ac)))=sub_1;
Matrix_1(Br(1):Br(size(Br)),Bc(1):Bc(size(Bc)))=sub_2;
Matrix_1(Cr(1):Cr(size(Cr)),Cc(1):Cc(size(Cc)))=sub_3;
Matrix_1(Dr(1):Dr(size(Dr)),Dc(1):Dc(size(Dc)))=sub_4;
Matrix_1(Er(1):Er(size(Er)),Ec(1):Ec(size(Ec)))=sub_5;
Matrix_1(Fr(1):Fr(size(Fr)),Fc(1):Fc(size(Fc)))=sub_6;
Matrix_1(Gr(1):Gr(size(Gr)),Gc(1):Gc(size(Gc)))=sub_7;
Matrix_1(Hr(1):Hr(size(Hr)),Hc(1):Hc(size(Hc)))=sub_8;
Matrix_1(Ir(1):Ir(size(Ir)),Ic(1):Ic(size(Ic)))=sub_9;
Matrix_1(Jr(1):Jr(size(Jr)),Jc(1):Jc(size(Jc)))=sub_10;
Matrix_1(Kr(1):Kr(size(Kr)),Kc(1):Kc(size(Kc)))=sub_11;
Matrix_1(Lr(1):Lr(size(Lr)),Lc(1):Lc(size(Lc)))=sub_12;
Matrix_1(Mr(1):Mr(size(Mr)),Mc(1):Mc(size(Mc)))=sub_13;
Matrix_1(Nr(1):Nr(size(Nr)),Nc(1):Nc(size(Nc)))=sub_14;
Matrix_1(Or(1):Or(size(Or)),Oc(1):Oc(size(Oc)))=sub_15;
Matrix_1(Pr(1):Pr(size(Pr)),Pc(1):Pc(size(Pc)))=sub_16;
end