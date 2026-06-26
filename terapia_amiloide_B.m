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
%supose que influix beta en tau lo mateix que al reves
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
%
xi_1=0.025; xi_2=0.025;
sigma_1=0.025; sigma_2=0.025;
%terapies
t_rx=20;%temps aplicacio terapia
T_vals=[115 250];   
Ti=0;
chi20=[0.005 0.01 0.025];
t_control=cell(length(T_vals),length(chi20));
Z_control=cell(length(T_vals),length(chi20));
t_terapia=cell(length(T_vals),length(chi20));
Z_terapia=cell(length(T_vals),length(chi20));

%% 3. SIMULACIONS
for iq=1:length(chi20)
chi2=chi20(iq);
for ir=1:length(T_vals)
    T=T_vals(ir);
    tspan=[Ti T];
    
    %control: 
    f_control=@(t,Z) model7(t,Z,omega_1,pi_0,...
        lambda_1,lambda_2,lambda_3,beta_0,...
        xi_1,xi_2,sigma_1,sigma_2,gamma_1,gamma_2,tau_0,...
        delta_1,delta_2,alpha_0,rho_1,rho_2,rho_3,nu_0);
    [t_control{ir,iq}, Z_control{ir,iq}]=ode45(f_control, tspan, Z_0);
    
    %terapia:
    f_terapia=@(t,Z) model7_terapia_amiloide_B(t,Z,omega_1,pi_0,...
        lambda_1,lambda_2,lambda_3,beta_0,...
        xi_1,xi_2,sigma_1,sigma_2,gamma_1,gamma_2,tau_0,...
        delta_1,delta_2,alpha_0,rho_1,rho_2,rho_3,nu_0,t_rx,chi2);
    [t_terapia{ir,iq}, Z_terapia{ir,iq}]=ode45(f_terapia, tspan, Z_0);
end

% 4. Gràfiques
noms_variables={'P', 'B', 'Q','C', 'T', 'A', 'N'};
for it=1:length(T_vals)
    T=T_vals(it);
    
    figure
    colororder(lines(7)) 
    hold on
    grid on
    
    %control
    plot(t_control{it,iq}, Z_control{it,iq}(:,[1 2 5 6 7]), 'LineWidth', 1.2, 'LineStyle', '--', 'HandleVisibility', 'off') 
    set(gca, 'ColorOrderIndex', 1)
    
    %terapia
    plot(t_terapia{it,iq}, Z_terapia{it,iq}(:,[1 2 5 6 7]), 'LineWidth', 1.8) 
    xline(t_rx, '-', 'LineWidth', 1.5);
    legend({'$P$','$B$','$T$','$A$','$N$'}, ...
           'Interpreter', 'latex', ...
           'Location', 'northwest')       
    xlim([Ti T])
    ylim([0 1.1])
    xlabel('Temps (anys)')
    ylabel('Biomarcadors')
    title(sprintf('Simulació teràpia anti-oligòmers (\\chi_2=%g)',chi2))
end
end
%% 4. Gràfiques
for k=[1 2 5 6 7]
for it=1:length(T_vals)
    T=T_vals(it);    
    figure
    colororder(lines(7)) 
    hold on
    grid on
    %control
    plot(t_control{it,1}, Z_control{it,1}(:,k), 'LineWidth', 1.2,...
        'LineStyle', '--', 'DisplayName', 'Control') ;
for iq=1:length(chi20)    
    %terapia
    plot(t_terapia{it,iq}, Z_terapia{it,iq}(:,k), 'LineWidth', 1.8,...
        'DisplayName', sprintf('\\chi_2=%g',chi20(iq))) ;      
end
xline(t_rx, '-', 'LineWidth', 1.5,'HandleVisibility','off');
xlim([Ti T])
ylim([0 1.1])
xlabel('Temps (anys)')
ylabel('Biomarcadors')
title(sprintf('Efecte de terapia anti-oligòmers sobre: %s', noms_variables{k}))
legend('Location', 'northwest')
end
end