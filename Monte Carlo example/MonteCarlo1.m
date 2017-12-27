%% Muestreando distribuciones conocidas
% extraído de Winkelried (2016)
% By Alex Carrasco, 2017

clear; clc; format bank; profile clear;
profile on;
%% I) set up
% Parametros
maxL=5;
L=2:maxL;
% Parametros
a=1; b=2;
mu = 2; sigma2 =1;
tdf = 3;
x2df = 6;
lambda=5;
alpha=3; beta=2;
qq = .9; %quantile de interés


% Definicion de distribuciones conocidas
dist = {'Unif'; 'Norm'; 't'; 'Chi2'; 'Exp'; 'Beta'};
tab1 = {'media';'var.';['Q' num2str(qq*100)];'4.EEN'};
tab2 = {['U(' num2str(a) ',' num2str(b) ')']; ['N(' num2str(mu) ',' num2str(sigma2) ')'] ; ...
        ['t(' num2str(tdf) ')']; ['X2(' num2str(x2df) ')']; ['exp(' num2str(lambda) ')'];...
        ['Beta(' num2str(alpha) ',' num2str(beta) ')']};
        
fda  = {@(x) a+(b-a)*rand(x,1); @(x) mu+sqrt(sigma2)*randn(x,1);  ...
         @(x) trnd(tdf,x,1); @(x) chi2rnd(x2df,x,1); @(x) exprnd(lambda,x,1); ...
         @(x) betarnd(alpha,beta,x,1)}; 
ppob   = [[(a+b)/2; (b-a)^2/12; qq], [mu;sigma2;norminv(qq,mu,sqrt(sigma2))],...
            [0;tdf/(tdf-2);icdf('t',qq,tdf)],[x2df;2*x2df;icdf('chi2',qq,tdf)], ...
            [lambda;lambda^2;icdf('exp',qq,lambda)], ...
            [alpha/(alpha+beta);alpha*beta/((alpha+beta)^2*(alpha+beta+1));icdf('beta',qq,alpha,beta)]];

%% II) Simulacion
C={};
for i=1:numel(dist)
    R.(dist{i})=[];
    Z=fda{i}(10^maxL);
    for j=L
       z=Z(1:10^j);
       R.(dist{i})=[R.(dist{i}),[mean(z);var(z,1);quantile(z,.9);4*std(z,1)/sqrt(10^j)]];   
    end
    aux1 = [dist(i);' ';' ';' ']; aux2=[tab2(i);' ';' ';' '];   % Con fines de reportar
    C = [C;[aux1,aux2,tab1,num2cell(R.(dist{i})), [num2cell(ppob(:,i));nan]]];     % Con fines de reportar
end

%% III) Reportando
T = cell2table(C,'VariableNames',{'Distrib','cont','cont_','M2','M3','M4','M5','Pob'});
T
profile viewer






