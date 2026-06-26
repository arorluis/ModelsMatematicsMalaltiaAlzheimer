function dZ_dt=modelmaji(t,Z,alpha,mu,lambda_0,lambda_1,rho,lambda_3, ...
    gamma_0,gamma_1,gamma_2,delta_0,delta_1,eta_0,eta_1,eta_2,n)
%VARIABLES
P=Z(1);
B=Z(2);
T=Z(3);
A=Z(4);
N=Z(5);

% EQUACIONS
dP_dt=alpha*B^n-mu*P;
dB_dt=lambda_0*B*(1-lambda_1*B)-n*alpha*B^n-rho*B*P+lambda_3*A;
dT_dt=gamma_0*T*(1-gamma_1*T)+gamma_2*B*T;
dA_dt=delta_0*B-delta_1*A;
dN_dt=(eta_0*T+eta_1*A)*(1-eta_2*N);

dZ_dt = [dP_dt; dB_dt; dT_dt; dA_dt; dN_dt];
end