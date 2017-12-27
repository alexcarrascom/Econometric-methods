%% MonteCarlo 5, regresión en series estacionarias
% extraído de Winkelried (2016)
% By Alex Carrasco, 2017
tic
clear; clc, close all;
%% [I] Inicializando
M = 10^4;                                       % n° de simulaciones [10^6]
Tsize = [50, 250, 350, 500]; maxT=max(Tsize);   % tamaño de muestra de analisis
T0 = 50;                                        % Muestra inicial para aproximar 'pasado remoto'

% Coeficientes
phi1 = 0.75;
phi2 = -0.2;
mu = 0;
df=3;

y = nan(maxT+T0,2); y(1:2,:)=mu;
phi  = nan(M,2,numel(Tsize));   % Media
t_p = nan(M,3,numel(Tsize));   % Media


%% [II] Comienza simulacion
for m=1:M
    e   = trnd(df,maxT+T0,1)/sqrt(df/(df-2));  % Chisquare con varianza normalizada;
    for t=3:maxT+T0
        y(t,1) = (1-phi1)*mu+phi1*y(t-1,1)+e(t);    % Proceso MA(1)
        y(t,2)=(1-phi1-phi2)*mu+phi1*y(t-1,2)+phi2*y(t-2,2)+e(t);    % Proceso AR(1)
    end
    y = y(T0-1:end,:);
    e = e(T0-1:end,:);
    
    for i=1:numel(Tsize)
        T=Tsize(i);
        Y = y(1:T+2,:);
        U = e(3:T+2);
        
        Y1 = Y(3:end,1); X1 = Y(2:end-1,1);
        Y2 = Y(3:end,2); X2 = [Y(2:end-1,2) Y(1:end-2,2)]; 
        % Proc 1
        phi_1=(X1'*X1)\(X1'*Y1);  s1 = T^(-1)*norm(Y1-phi_1*X1)^2;  
        % Proc 2
        aux1=inv(X2'*X2);
        aux2=(X2'*X2)\(X2'*Y2); s2 = T^(-1)*norm(Y2-X2*aux2)^2;
        
        phi(m,1,i) = phi_1;
        phi(m,2,i) = aux2(1);
        t_p(m,1,i) = (phi(m,1,i)-phi1)/(sqrt(s1/(X1'*X1)));
        t_p(m,2,i) = (phi(m,2,i)-phi1)/(sqrt(s2*aux1(1)));
    end
end

%% [III] Graficos
styles={'-.b',':b','--b','b'};
xn = (-4:0.01:4)';
yn = (2*pi)^(-1/2)*exp(-.5*xn.^2);
titul={'AR(1)','AR(2)'};

figure(1)
for k=1:2 % columnas
    if k==1
        for i=0:1
            subplot(2,2,i*2+k)
            for j=1:numel(Tsize)
                [y,~,x]=mykernel(phi(:,i+1,j),'kernel',2,'eachseries');
                plot(x,y,styles{j});
                hold on
            end
            xlim([0.4 1]);
            title([titul{i+1} ', \phi_1']);
        end
        legend(['T=' num2str(Tsize(1))],['T=' num2str(Tsize(2))], ...
            ['T=' num2str(Tsize(3))],['T=' num2str(Tsize(4))],'location','best');
    else 
        for i=0:1
            subplot(2,2,i*2+k)
            for j=1:2
                [y,~,x]=mykernel(t_p(:,i+1,j),'kernel',2,'eachseries');
                plot(x,y,styles{j*2});
                hold on
            end
            plot(xn,yn,'-b','linewidth',1.5);
            xlim([-3 3]);
            ylim([0 0.48]);
            title([titul{i+1} ', t_{\phi}']);
        end
        legend(['T=' num2str(Tsize(1))],['T=' num2str(Tsize(2))],'N(0,1)','location','best');
    end
end

% print -depsc2 MonteCarlo5.eps;
% eps2pdf('MonteCarlo5.eps');
 
delete *.eps;
toc










