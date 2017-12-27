%% Montecarlo 4: Leyes asintóticas, variables con dependencia
% extraído de Winkelried (2016)
% By Alex Carrasco, 2017
tic
clear all; clc; close all;

%% [I] Inicializando
M = 10^3;                                       % n° de simulaciones [10^5]
Tsize = [50, 250, 350, 500]; maxT=max(Tsize);   % tamaño de muestra de analisis
T0 = 50;                                        % Muestra inicial para aproximar 'pasado remoto'
% Coeficientes de los procesos
phi   = 0.75;
theta = 0.5;
mu = 1;
df =3;

y = nan(maxT+T0,3); y(1,:)=mu;
ybar = nan(M,3,numel(Tsize));   % Media
stdy = nan(M,3,numel(Tsize));   % Media
corr1 = nan(M,3,numel(Tsize));   % Coeficiente de correlacion

%% [II] Comienza simulacion
for m=1:M
    e   = trnd(df,maxT+T0,1)/sqrt(df/(df-2));  % Chisquare con varianza normalizada;
    for t=2:maxT+T0
        y(t,1) = mu+e(t)+theta*e(t-1);    % Proceso MA(1)
        y(t,2)=(1-phi)*mu+phi*y(t-1,2)+e(t);    % Proceso AR(1)
        y(t,3)=(1-phi)*mu+phi*y(t-1,3)+e(t) + theta*e(t-1); % ARMA(1,1)
    end
    y = y(T0+1:end,:);
    
    for i=1:numel(Tsize)
        T=Tsize(i);
        Y = y(1:T,:);
        ybar(m,:,i) = mean(Y);
        stdy(m,:,i) = std(Y);
        corr1(m,:,i) = (T-1)^(-1)*diag( ( Y(2:end,:)-kron(ones(T-1,1),mean(Y)) )'* ...
                        ( Y(1:end-1,:)-kron(ones(T-1,1),mean(Y)) ) )'./var(Y);
    end
end

D = diag([1/(1+theta),1-phi,(1-phi)/(1+theta)]);
zbar = [];
zbar(:,:,1) = sqrt(Tsize(1))*(ybar(:,:,1)-mu)*D; 
zbar(:,:,2) = sqrt(Tsize(2))*(ybar(:,:,2)-mu)*D; 
zbar(:,:,3) = sqrt(Tsize(3))*(ybar(:,:,3)-mu)*D; 
zbar(:,:,4) = sqrt(Tsize(4))*(ybar(:,:,4)-mu)*D; 

%% [III] Graficos
styles={'-.b',':b','--b','b'};
xn = (-4:0.01:4)';
yn = (2*pi)^(-1/2)*exp(-.5*xn.^2);
titul={'MA(1)','AR(1)','ARMA(1,1)'};

figure(1)
for k=1:3
    if k==1
        for i=0:2
            subplot(3,3,i*3+k)
            for j=1:numel(Tsize)
                [y,~,x]=mykernel(ybar(:,i+1,j),'kernel',2,'eachseries');
                plot(x,y,styles{j});
                hold on
            end
            xlim([mu-2 mu+2]);
            title([titul{i+1} ', ybar']);
        end
%         legend(['T=' num2str(Tsize(1))],['T=' num2str(Tsize(2))], ...
%             ['T=' num2str(Tsize(3))],['T=' num2str(Tsize(4))],'location','best');
    elseif k==2
        for i=0:2
            subplot(3,3,i*3+k)
            for j=1:2
                [y,~,x]=mykernel(zbar(:,i+1,j),'kernel',2,'eachseries');
                plot(x,y,styles{j*2});
                hold on
            end
            plot(xn,yn,'-b','linewidth',1.5);
            xlim([-4 4]);
            ylim([0 0.5]);
            title([titul{i+1} ', zbar']);
        end
        legend(['T=' num2str(Tsize(1))],['T=' num2str(Tsize(2))],'location','best');
    else
        for i=0:2
            subplot(3,3,i*3+k)
            for j=1:numel(Tsize)
                [y,~,x]=mykernel(corr1(:,i+1,j),'kernel',2,'eachseries');
                plot(x,y,styles{j});
                hold on
            end
            xlim([0 1]);
            title([titul{i+1} ', \rho_1']);
        end
        legend(['T=' num2str(Tsize(1))],['T=' num2str(Tsize(2))], ...
            ['T=' num2str(Tsize(3))],['T=' num2str(Tsize(4))],'location','best');        
    end
end

% print -depsc2 MonteCarlo4.eps;
% eps2pdf('MonteCarlo4.eps');

delete *.eps;

toc