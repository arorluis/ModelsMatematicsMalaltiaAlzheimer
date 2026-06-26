clear; clc; close all;
% Paràmetres
xi_1=0.025;
xi_2=0.025;
sigma_1=0.025;
sigma_2=0.025;
% Condicions inicials
Q0=0.5;
C0=0;
t=linspace(0,115,1000);
% Solució cas igualtat
Q=Q0.*exp(-xi_2*t)+(xi_1/xi_2).*(1-exp(-xi_2*t));
C=C0.*exp(-sigma_2*t)+ sigma_1*(Q0-xi_1/xi_2).*t.*exp(-sigma_2*t) ...
   +(sigma_1/sigma_2)*(xi_1/xi_2).*(1-exp(-sigma_2*t));
% Gráfica
figure
plot(t,Q,'LineWidth',2)
hold on
plot(t,C,'LineWidth',2)
xlabel('temps)')
ylabel('Biomarcadors')
legend('Q(t)','C(t)','Location','best')
grid on
title(sprintf('Solució analítica cas \\sigma_2=%g, \\xi_2=%g',sigma_2, xi_2))
xlim([0 115])
%% Cas desigualtat
clear; clc; 
% Paràmetres
xi_1=0.025;
xi_2=0.025;
sigma_1=0.025;
sigma_2=0.03;
% Condicions inicials
Q0=0.5;
C0=0;
t=linspace(0,115,1000);
% Solució cas igualtat
Q=Q0.*exp(-xi_2*t)+(xi_1/xi_2).*(1-exp(-xi_2*t));
C=C0.*exp(-sigma_2.*t) +...
    (sigma_1/(sigma_2-xi_2)).*(Q0-xi_1/xi_2).*( exp(-xi_2*t)-exp(-sigma_2*t))+...
    (sigma_1/sigma_2)*(xi_1/xi_2).*(1-exp(-sigma_2*t));
% Gráfica
figure
plot(t,Q,'LineWidth',2)
hold on
plot(t,C,'LineWidth',2)
xlabel('temps')
ylabel('Biomarcadors')
legend('Q(t)','C(t)','Location','best')
grid on
title(sprintf('Solució analítica cas \\sigma_2=%g, \\xi_2=%g',sigma_2, xi_2))
xlim([0 115])