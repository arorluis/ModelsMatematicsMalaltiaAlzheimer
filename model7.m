function dZ_dt = model7(t,Z,omega_1,pi_0,...
    lambda_1,lambda_2,lambda_3,beta_0,...
    xi_1,xi_2,sigma_1,sigma_2,gamma_1,gamma_2,tau_0,...
    delta_1,delta_2,alpha_0,rho_1,rho_2,rho_3,nu_0)
% VARIABLES
P = Z(1);
B = Z(2);
Q = Z(3);
C = Z(4);
T = Z(5);
A = Z(6);
N = Z(7);

% EQUACIONS
dP_dt = omega_1*B*(1 - pi_0*P);
dB_dt = (lambda_1*A + lambda_2*C - lambda_3*P)*(1 - beta_0*B);
dQ_dt = xi_1 - xi_2*Q;
dC_dt = sigma_1*Q - sigma_2*C;
dT_dt = (gamma_1*B + gamma_2*A)*(1 - tau_0*T);
dA_dt = (delta_1*B + delta_2*T)*(1 - alpha_0*A);
dN_dt = (rho_1*T + rho_2*A + rho_3*P) *(1 - nu_0*N);

dZ_dt = [dP_dt;dB_dt;dQ_dt;dC_dt;dT_dt;dA_dt;dN_dt];
end