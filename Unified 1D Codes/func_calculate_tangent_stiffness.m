function [J_el] = func_calculate_tangent_stiffness(Damage,dDdestar,B,u_current_element,ConstitutiveMatrix,el,weight_integrationpoint,n_IntegrationPoints,J)
% This function calculates the tangent stiffness matrix of the element
% under consideration

% Initialize tangent stiffness matrix components
[J11,J12,J21,J22]=deal(0);


% Calculating the Element tangent stiffness matrix
% J11=dR1/du1

A1=(-1*(dDdestar)*B(1)*(B(1)*ConstitutiveMatrix(1,el)*B(1))*u_current_element(1));
A2=(1-Damage)*B(1)*ConstitutiveMatrix(1,el)*B(1);
A3=(-1*(dDdestar)*B(1)*(B(1)*ConstitutiveMatrix(1,el)*B(2))*u_current_element(2));

J11=J11+((A1+A2+A3)*weight_integrationpoint(n_IntegrationPoints)*J);

% J12=dR1/du2

B1=(-1*(dDdestar)*B(2)*(B(1)*ConstitutiveMatrix(1,el)*B(1))*u_current_element(1));
B2=(1-Damage)*B(1)*ConstitutiveMatrix(1,el)*B(2);
B3=(-1*(dDdestar)*B(2)*(B(1)*ConstitutiveMatrix(1,el)*B(2))*u_current_element(2));

J12=J12+((B1+B2+B3)*weight_integrationpoint(n_IntegrationPoints)*J);

% J21=dR2/du1

C1=(-1*(dDdestar)*B(1)*(B(2)*ConstitutiveMatrix(1,el)*B(1))*u_current_element(1));
C2=(1-Damage)*B(2)*ConstitutiveMatrix(1,el)*B(1);
C3=(-1*(dDdestar)*B(1)*(B(2)*ConstitutiveMatrix(1,el)*B(2))*u_current_element(2));

J21=J21+((C1+C2+C3)*weight_integrationpoint(n_IntegrationPoints)*J);

% J22=dR2/du2

D1=(-1*(dDdestar)*B(2)*(B(2)*ConstitutiveMatrix(1,el)*B(1))*u_current_element(1));
D2=(1-Damage)*B(2)*ConstitutiveMatrix(1,el)*B(2);
D3=(-1*(dDdestar)*B(2)*(B(2)*ConstitutiveMatrix(1,el)*B(2))*u_current_element(2));

J22=J22+((D1+D2+D3)*weight_integrationpoint(n_IntegrationPoints)*J);

J_el=[J11,J12;J21,J22];

end