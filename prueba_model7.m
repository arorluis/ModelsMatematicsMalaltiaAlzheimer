clear;
clc;
close all;

%%
%  VALORS INICIALS
P_0=0;
B_0=0.01;
Q_0=0.5;
C_0=0;%al principi no hi ha C99
T_0=0;
A_0=0;
N_0=0;
Z_0=[P_0; B_0; Q_0; C_0; T_0; A_0; N_0];

% PARÀMETRES
omega_1=0.09;
gamma_1=0.025;
gamma_2=0.025;
lambda_1=0.025;lambda_2=0.025;lambda_3=0.025;
pi_0=1;
beta_0=1;
tau_0=1;
alpha_0=1;
nu_0=1;
%inspirats petrella:
%supose que la mes influent es la tau
rho_1=0.025;
rho_2=0.010;
rho_3=0.010;
%la neuroinflamació s'activa abans i potencia la tau més que al revés
delta_1= 0.05;delta_2= 0.01;
%subsistema
xi_1=0.025; xi_2=0.025;
sigma_1_0=[0 0.025]; sigma_2=0.025;
%sense terapies
T_vals=[115 250];
Ti=0;
t=cell(1, length(T_vals));
Z=cell(1, length(T_vals));
for isigma=1:length(sigma_1_0)
sigma_1=sigma_1_0(isigma);
%% MODEL
f=@(t,Z) model7(t,Z,omega_1,pi_0,...
    lambda_1,lambda_2,lambda_3,beta_0,...
    xi_1,xi_2,sigma_1,sigma_2,gamma_1,gamma_2,tau_0,...
    delta_1,delta_2,alpha_0,rho_1,rho_2,rho_3,nu_0);
%% SIMULACIONS
for it=1:length(T_vals)
    T=T_vals(it);
    tspan=[Ti T];
    [t{it}, Z{it}]=ode45(f, tspan, Z_0);
end
%% GRÀFIQUES
for it=1:length(T_vals)
    figure
    colororder(lines(7))
    hold on
    grid on
    plot(t{it}, Z{it}, 'LineWidth', 1.5)
    legend({'$P$','$B$','$Q$','$C$','$T$','$A$','$N$'}, ...
           'Interpreter','latex','Location','best')
    xlabel('Temps (anys)')
    ylabel('Biomarcadors')
    xlim([Ti T_vals(it)])
    ylim([0,1.1]);
    title(sprintf('Evolució temporal (\\sigma_1=%g)',sigma_1))
end
end