%% Montecarlo 2: Eficiencia relativa de dos estimadores
% extraído de Winkelried (2016)
% By Alex Carrasco, 2017

% ER(theta1,theta2)=var(theta1)/var(theta2)
clear; clc; tic;
% Definamos el tamaño de muestra y el numero de replicaciones
n = 100;
M = 10^5;

% Considere una distribucion t-student con v grados de libertad
maxV = 20;
ER=nan(maxV,1);
for v=1:maxV
    t_v=trnd(v,n,M);
    ER(v)=var(mean(t_v),1)/var(median(t_v),1);
end
toc

%% Grafico
plot(ER,'o-','linewidth',1.5,'markersize',9); hold on;
plot([1 maxV],[1 1],'k','linewidth',1.25);
plot([1 maxV],[2/pi 2/pi],'k--','linewidth',1.25); 
axis tight; ylim([0 1.8]);
xlabel('\upsilon - grados de libertad'); ylabel('ER: media vs mediana')

% print -depsc2 MC2ER.eps;
% eps2pdf('MC2ER.eps'); 
delete *.eps;

