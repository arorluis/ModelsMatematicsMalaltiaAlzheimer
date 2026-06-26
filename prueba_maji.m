clear;
clc;
close all;

%%
%  VALORS INICIALS
n0=[1 2 3 4 5];
%n0=1;
for in=1:length(n0)
n=n0(in);
P_0=8.4;
B_0=2.6;
T_0=4.6;
A_0=6.8;
N_0=0.5;
Z_0=[P_0; B_0; T_0; A_0; N_0];

% PARÀMETRES
alpha=0.8;
mu=0.035;
lambda_0=2.1;
lambda_1=1.09;
rho=0.9;
lambda_3=1.26;
gamma_0=1.2;
gamma_1=0.8;
gamma_2=0.9;
delta_0=0.8;
delta_1=0.095;
eta_0=6.9e-3;
eta_1=2e-3;
eta_2=1;

T_vals=[150];
Ti=0;

t=cell(1, length(T_vals));
Z=cell(1, length(T_vals));

%% MODEL
f=@(t,Z) modelmaji(t, Z, alpha, mu, lambda_0, lambda_1, rho, lambda_3, ...
    gamma_0, gamma_1, gamma_2, delta_0, delta_1, eta_0, eta_1, eta_2, n);

%% SIMULACIONS
for it=1:length(T_vals)
    T=T_vals(it);
    tspan=[Ti T];
    [t{it}, Z{it}]=ode45(f, tspan, Z_0);
end

%% GRÀFIQUES

for it=1:length(T_vals)
figure
colororder(lines(5))
hold on
grid on
plot(t{it}, Z{it}, 'LineWidth', 1.5)
legend({'$P$','$B$','$T$','$A$','$N$'}, ...
       'Interpreter','latex', ...
       'Location','northeast')

xlabel('Temps (dies)')
ylabel('Biomarcadors')
xlim([Ti T_vals(it)])
title(sprintf('Evolució temporal n=%d',n))
end

%% Normalitzem per a vore millor
for it=1:length(T_vals)
figure
colororder(lines(5))
hold on
grid on
Znorm=Z{it} ./ max(Z{it});
plot(t{it}, Znorm, 'LineWidth',1.5)
legend({'$P$','$B$','$T$','$A$','$N$'}, ...
       'Interpreter','latex', ...
       'Location','best')

xlabel('Temps (dies)')
ylabel('Valor normalitzat')
title(sprintf('Variables normalitzades n=%d',n))
ylim([0,1.2]);
end

end