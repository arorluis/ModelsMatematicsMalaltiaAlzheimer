clear;
clc;
close all;

%% VALORS INICIALS
P_0=0;
B_0=0.01;
Q_0=0.5;
C_0=0; % al principi no hi ha C99
T_0=0;
A_0=0;
N_0=0;
Z_0=[P_0; B_0; Q_0; C_0; T_0; A_0; N_0];

%% PARÀMETRES
% del petrella
omega_1=0.09;
gamma_1=0.025;
gamma_2=0.025;

% capacitats
pi_0=1;
beta_0=1;
tau_0=1;
alpha_0=1;
nu_0=1;

% inspirats petrella
rho_1=0.025;
rho_2=0.010;
rho_3=0.010;

% la neuroinflamació s'activa abans i potencia la tau més que al revés
delta_1=0.05;
delta_2=0.01;

% subsistema
xi_1=0.025;
xi_2=0.025;
sigma_1=0.025;
sigma_2=0.025;

% 
lambda_1=0.025;
lambda_2=0.025;

% escombratge en lambda3 incloent el equilibri
lambda_3_eq=(pi_0*(alpha_0*lambda_2*sigma_1*xi_1 + ...
    lambda_1*sigma_2*xi_2))/(alpha_0*sigma_2*xi_2);

v1=linspace(0, lambda_3_eq, 5);
v2=linspace(lambda_3_eq, 2*lambda_3_eq, 5);

lambda_3_0=[v1, v2(2:end)];

t_terapia=1e6;

T_vals=[250 1e3];
Ti=0;

%% VARIABLES A REPRESENTAR
nom_vars={'P','B','Q','C','T','A','N'};

%% SIMULACIONS 
t=cell(1,length(lambda_3_0));
Z=cell(1,length(lambda_3_0));
for it=1:length(T_vals)
T=T_vals(it);
tspan=[Ti T];

for ilambda=1:length(lambda_3_0)

    lambda_3=lambda_3_0(ilambda);

    f=@(t,Z) model7(t,Z,omega_1,pi_0,...
        lambda_1,lambda_2,lambda_3,beta_0,...
        xi_1,xi_2,sigma_1,sigma_2,gamma_1,gamma_2,tau_0,...
        delta_1,delta_2,alpha_0,rho_1,rho_2,rho_3,nu_0);

    [t{it,ilambda}, Z{it,ilambda}]=ode45(f,tspan,Z_0);

end
end
%% GRÀFIQUES
for it=1:length(T_vals)
T=T_vals(it);
for ivar=1:7

    figure
    hold on
    grid on

    for ilambda=1:length(lambda_3_0)

        plot(t{it,ilambda}, ...
             Z{it,ilambda}(:,ivar), ...
             'LineWidth',1.5,...
             'DisplayName',sprintf('$\\lambda_3/\\bar{\\lambda}_3=%.2f$',...
             lambda_3_0(ilambda)/lambda_3_eq));

    end

    xlabel('Temps (anys)')
    ylabel('Biomarcadors')
    xlim([Ti T])
    ylim([-0.2 1.1])

    title(sprintf('Canvi de $%s(t)$ amb $\\lambda_3$', ...
        nom_vars{ivar}), ...
        'Interpreter','latex')

    legend('show','Interpreter','latex','Location','best')

end
end