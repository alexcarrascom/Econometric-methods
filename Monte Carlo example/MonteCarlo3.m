%% Monte Carlo 3: Estimador hibrido (Pre testing)
% extraído de Winkelried (2016)
% By Alex Carrasco, 2017

% Regresores no estocásticos
clear all; clc; format long; profile clear;
rng('default');
%% [I] configuracion de la simulacion
tic;
n=20; 
M = 10^5;
k = 1;  % k regresores
b=1;
X_nc = 5+3*randn(n,k);
X  = [ones(n,1),X_nc];  
iX2_nc = inv(X_nc'*X_nc);
iX2 = inv(X'*X);
sigma = 2;

%% [II] Configuracion de los distintos DGP
alphamin=0; alphamax=3; nalpha=20;
alphas=linspace(alphamin,alphamax,nalpha);
t_crit = icdf('t',.975,n-k-1);
beta_ols=nan(k+1,M,nalpha);
beta_mcr=nan(k,M,nalpha);
beta_pt=nan(k,M,nalpha);
rechazo = nan(M,nalpha);

%% [III] simulacion
for jj=1:numel(alphas)
    alpha = alphas(jj);
    beta0  =[alpha; b*ones(k,1)];

    for m=1:M
        y = X*beta0 + sigma*randn(n,1);
        beta_ols(:,m,jj) =  iX2*X'*y; sigma2est = norm(y-X*beta_ols(:,m,jj))^2/(n-k-1);
        t_a   = beta_ols(1,m,jj)/(sigma2est*iX2(1,1));
        beta_mcr(:,m,jj) = iX2_nc*X_nc'*y;
        rechazo(m,jj) = (abs(t_a)>t_crit);
        beta_pt(:,m,jj)  = beta_mcr(:,m,jj)- rechazo(m,jj)*(beta_mcr(:,m,jj)-beta_ols(2:k+1,m,jj));
    end
end

%% [IV] GRAFICO
alpplot=3:3:18;

% a) Prob
betaols=squeeze(beta_ols(2,:,:));
betamcr=squeeze(beta_mcr);
betapt=squeeze(beta_pt);
B(:,:,1) = betaols; B(:,:,2) = betamcr; B(:,:,3) = betapt;
poderf = mean(rechazo);
media =squeeze(mean(B)); 
sesgob2=media-beta0(2);
varb2  = squeeze(var(B,1));
ECMb2  = varb2+ sesgob2.^2;

figure(1);
subplot(2,2,1)
plot(alphas,poderf,'linewidth',1.5,'linestyle','--'); axis tight; title('Función poder');
subplot(2,2,2)
h=plot(alphas,sesgob2); axis tight; ylim([-0.01,0.2]); title('Sesgo');
set(h(1),'color','b'); set(h(2),'color','b','linestyle','--');
set(h(3),'color','b','linewidth',2);
legend('MCO','MCR','PT');
subplot(2,2,3)
h=plot(alphas,varb2); axis tight; ylim([0,0.03]); title('varianza');
set(h(1),'color','b'); set(h(2),'color','b','linestyle','--');
set(h(3),'color','b','linewidth',2);
legend('MCO','MCR','PT','location','northwest');
subplot(2,2,4)
h=plot(alphas,ECMb2); axis tight; ylim([-0.01,0.1]); title('ECM');
set(h(1),'color','b'); set(h(2),'color','b','linestyle','--');
set(h(3),'color','b','linewidth',2);
legend('MCO','MCR','PT','location','northwest');

% print -depsc2 'MC21.eps';
% eps2pdf('MC21.eps');

%% b)
cc=1;
figure(2); 
for kk=alpplot
    beta2ols=betaols(:,kk);
    beta2mcr=betamcr(:,kk);
    beta2pt=betapt(:,kk);
    [yols,~,Xols]=mykernel(beta2ols,'kernel',2,'eachseries');
    [ymcr,~,Xmcr]=mykernel(beta2mcr,'kernel',2,'eachseries');
    [ypt,~,Xpt]=mykernel(beta2pt,'kernel',2,'eachseries');
    
    subplot(3,2,cc); hold on;
    plot(Xols,yols,'-');
    plot(Xmcr,ymcr,':','linewidth',1.25);
    plot(Xpt,ypt,'-','linewidth',2);
    legend('MCO','MCR','PT','location','northwest');
    title(['\alpha=' num2str(alphas(kk),2) ', poder=' num2str(poderf(kk),2)]);  
    cc=1+cc;    
end

% print -depsc2 'MC22.eps';
% eps2pdf('MC22.eps');

% end
toc
delete *.eps;







