%% MonteCarlo 6: Test de raíz unitaria
% extraído de Winkelried (2016)
% By Alex Carrasco, 2017
tic
clear all; clc; close all;
rng(1115);
%% [I] Inicializando
M = 10^4; % n° simulaciones (en el documento se usó 10^6)
T = 10^2;
T0 = 40;

alpha=0;
df=3;

y0=0;
rho  = nan(M,3);  
t_T = nan(M,3);   

%% [II] Comienza simulacion
for m=1:M
    e   = trnd(df,T+T0,1)/sqrt(df/(df-2));  % t-student con varianza normalizada;
    cse = cumsum(e);
    y = cse + y0;                          
    Y = y(T0+1:end);        % datos observados   
    X1 = y(T0:end-1);
    X2 = [ones(T,1) X1];
    X3 = [X2 cumsum(ones(T,1))];
    
    % Caso 1
    rho(m,1) = X1'*Y/(X1'*X1); s2 = T^(-1)*norm(Y-X1*rho(m,1))^2;
    t_T(m,1) = (rho(m,1)-1)/sqrt(s2/(X1'*X1));
    
    % Caso 2
    beta = (X2'*X2)\(X2'*Y); s2 = T^(-1)*norm(Y-X2*beta)^2;
    varbeta = s2*inv(X2'*X2);
    rho(m,2) = beta(2);
    t_T(m,2) = (beta(2)-1)/sqrt(varbeta(2,2));
    
    % Caso 2
    beta = (X3'*X3)\(X3'*Y); s2 = T^(-1)*norm(Y-X3*beta)^2;
    varbeta = s2*inv(X3'*X3);
    rho(m,3) = beta(2);
    t_T(m,3) = (beta(2)-1)/sqrt(varbeta(2,2));    
end

%% [III] Graficos
styles={'-b',':b','--b'};
xt = (-3:0.01:3)';

figure(1)
subplot(2,1,1)
for i=1:3
    [y,~,x]=mykernel(rho(:,i),'kernel',2,'eachseries');
    plot(x,y,styles{i},'linewidth',1.25);
    hold on;
end
xlim([.8 1.1]);
legend('Caso 1', 'Caso 2', 'Caso 3','location','best');

subplot(2,1,2)
for i=1:3
    [y,~,x]=mykernel(t_T(:,i),'kernel',2,'eachseries');
    plot(x,y,styles{i},'linewidth',1.25);
    hold on;
end
plot(xt,tpdf(xt,T),'b','linewidth',2); xlim([-4, 3]);
legend('Caso 1', 'Caso 2', 'Caso 3','t-student(T)','location','best');

% print -depsc2 MonteCarlo6.eps;
% eps2pdf('MonteCarlo6.eps');
 
delete *.eps;
toc

