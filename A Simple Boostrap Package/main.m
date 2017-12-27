%% BOOTSTRAP package for matlab
% Main driver to replicate examples
% Course: Econometrics II
% Prof.: Pedro CL Souza
% Student: Alex Carrasco
% ********************************************************
%   By Alex Carrasco, december 2017 (PUC Rio)  
% ********************************************************
tic;
clear all; clc; close all;
cc=1;

path_g = [pwd '\Graphs']; 
% imprime = @(x) print( gcf, '-depsc2', [path_g filesep x]);
% imprpdf = @(x) eps2pdf( [path_g filesep x '.eps']);
formataxis  = @(x) set(x, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 23, 'Box', 'Off', 'PlotBoxAspectRatio', [1 0.75 1]);
formatlegend  = @(x) set(legend, 'Location', x, 'Orientation', 'Vertical', 'Box', 'Off', 'Fontsize', 23, 'Fontangle', 'normal');
label_x   =@(x) xlabel(x,'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 20,'interpreter','latex');
label_y   =@(x) ylabel(x,'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 20,'interpreter','latex');

%% [I] Naive bootstrap
rng(1115);      % controlling random numbers generator (for reproducibility)
N     = 150;
B     = 3000;
grids = 100;

lambda=10;
T=@(x) [mean(x) var(x) var(x,1)];  % estimators

X=random('poisson',lambda,N,1);

% Printing results
[t,Y,EDF,Xgrid]=naive_boots(X,T,'n_replic',B,'n_grid',grids);
fprintf('\n======================================\n');
fprintf(' [I] Naive Bootstrap: estimating bias\n');
fprintf('======================================\n');
lambest = T(X)
bias    = mean(t)- T(X);
badj    = T(X)- bias

% Plot
figure(cc); hold on;
F = poisscdf(1:Xgrid(end),lambda);
h1=stairs(1:Xgrid(end),F);
h2=plot(Xgrid,EDF);
set(h1,'color',[0.7 0.7 0.7],'linestyle','-','linewidth',2.5);
set(h2,'color',[0 0 0.5],'linestyle','none','linewidth',2,'marker','o','markersize',4);
axis tight; xlim([0 Xgrid(end)])
set(gca,'xtick',0:2:Xgrid(end));
formataxis(gca);
legend('F_0','EDF');
formatlegend('Best');
% imprime('Ex1');
% imprpdf('Ex1'); 
cc=cc+1;

%% [II] Kernel densities
rng(1115);           % controlling random numbers generator (for reproducibility)
N  = 300;
S  = [1 .5; .5 1];   % covariance matrix
mu = 5*ones(2,1);    % mean
D  = mvnrnd(mu,S,N); % generating sample

% Uniform kernel estimator (default) 
[fu,Fu,Xgrid] = mykernel(D,'n_grids',100);

% Gaussian kernel estimator (kernel=2) 
[fn,Fn]  =  mykernel(D,'n_grids',100,'kernel',2);

% Epanechnikov kernel estimator (kernel=3) 
[fe,Fe]  =  mykernel(D,'n_grids',100,'kernel',3,'bandwith',1);

% Bivariate Gaussian PDF 
Si=inv(S);
mvrnd=@(x,y) (2*pi)^(-1)*(det(S))^(-1/2)*exp(-0.5*(Si(1,1)*(x-mu(1)).^2+Si(2,2)*(y-mu(2)).^2 ...
                                           +2*Si(1,2)*(x-mu(1)).*(y-mu(2)) )); 

% Plots
ww =  1;
dist = {'u','n','e'};
[X1,X2]=meshgrid(Xgrid(:,1),Xgrid(:,2));
ftrue= mvrnd(X1,X2);

figure(cc)
h=surf(X1,X2,ftrue); axis tight;
set(h,'edgealpha',.5,'facelighting','flat','facealpha',.75);
view(-35,30);
label_x('var1'); label_y('var2');
formataxis(gca);
% imprime(['Ex2_' num2str(ww)]);
% imprpdf(['Ex2_' num2str(ww)]);
cc=cc+1; ww=ww+1;

figure(cc); hold on;
h1 = plot(D(:,1),D(:,2));
h2 = contour(X1,X2,ftrue,15);
set(h1,'markersize',8,'marker','x','color','blue','linestyle','none');
label_x('var1'); label_y('var2');
formataxis(gca);
% imprime(['Ex2_' num2str(ww)]);
% imprpdf(['Ex2_' num2str(ww)]);
cc=cc+1; ww=ww+1;

for i=1:numel(dist)
    figure(cc)
    eval(['h=surf(X1,X2, transpose(f' dist{i} '));']); axis tight;
    set(h,'edgealpha',.5,'facelighting','flat','facealpha',.75);
    view(-35,30);
    label_x('var1'); label_y('var2');
    formataxis(gca); 
%     imprime(['Ex2_' num2str(ww)]);
%     imprpdf(['Ex2_' num2str(ww)]);
    cc=cc+1; ww=ww+1;
    
    figure(cc); hold on;
    h1 = plot(D(:,1),D(:,2));
    eval(['contour(X1,X2,transpose(f' dist{i} '),15);']);
    set(h1,'markersize',8,'marker','x','color','blue','linestyle','none');
    label_x('var1'); label_y('var2');    
    formataxis(gca);
    legend('Data','Isolines'); formatlegend('best');
%     imprime(['Ex2_' num2str(ww)]);
%     imprpdf(['Ex2_' num2str(ww)]);
    cc=cc+1; ww=ww+1;
    
    figure(cc);
    eval(['h=surf(X1,X2, transpose(F' dist{i} '));']); axis tight;
    set(h,'edgealpha',.5,'facelighting','flat','facealpha',.75);
    view(-35,30);
    label_x('var1'); label_y('var2');
    formataxis(gca); 
%     imprime(['Ex2_' num2str(ww)]);
%     imprpdf(['Ex2_' num2str(ww)]);
    cc=cc+1;  ww=ww+1;   
end

%% [III] Non-parametric bootstrap
load BOV_data.mat;       % IBOVESPA series: Dsp (series), ran (sample range)
p = [.05 .3 .5 .8 .95];  % percentiles of interest
T = @(x) [mean(x) quantile(x,p)];  

% smoothed bootstrap using Gaussian kernel
[t,Y,F,Xgrid,f]=nonparam_boots(Dsp,T,'kernel',2,'n_replic',1e3);

% Printing results
fprintf('\n============================================================\n');
fprintf(' [II] Non-parametric Bootstrap: estimating percentiles (IBOVESPA)\n');
fprintf('============================================================\n');
fprintf('Estimation using sample\n');
T(Dsp)
fprintf('Estimation using bootstrap sample\n');
mean(t)
fprintf('Estimation using bootstrap sample (adjusted)\n');
2*T(Dsp)-mean(t)

% plots
figure(cc)
h=plot(ran,Dsp); axis tight; hold on;
plot([ran(1) ran(end)],[0 0],':k')
set(h,'color',[0 0 .5],'linewidth',1.25);
set(gca,'xtick',ran(1):47:ran(end),'xticklabel',dat2str(ran(1):47:ran(end)));
formataxis(gca); cc=cc+1;
% imprime('Ex3_1');
% imprpdf('Ex3_1');

figure(cc)
h=plot(Xgrid,[f normpdf(Xgrid,0,6)]);
set(h(1),'color',[0 0 .5],'linewidth',2);
set(h(2),'color',[.5 0 0],'linewidth',2,'linestyle','none','marker','x');
axis tight;
formataxis(gca);
legend('Normal(0,6)','KDE - IBOVESPA');
formatlegend('northwest'); cc=cc+1;
% imprime('Ex3_2');
% imprpdf('Ex3_2');

%% [IV] Residual bootstrap (Panel data)
rng(1115);    % controlling random numbers generator (for reproducibility)  
N = 500;      
B = 200;
lambda = [0.1  0.5];    % parameters for each covariate
Sigma  = 15*[1 .5; .5 1]; % residuals variance matrix
btrue  = [5 10]';      % parameters of interest 

% generating random sample
x1 = exprnd(lambda(1),N,2);   % t=1 
x2 = exprnd(lambda(2),N,2);   % t=2
v  = mvnrnd([0; 0],Sigma,N);

% observable sample
X  = reshape([x1'; x2'],2,N*2)';
y  = X*btrue+vec(v');

% POLS
[t_p,St_p,V_p,bhat_p,S_p,D_p]=residual_boots(y,X,N,'n_replic',B);
Sbeta_p  = (X'*X)\((X'*kron(eye(N),Sigma)*X)/(X'*X));
aux      = bsxfun(@minus,t_p,mean(t_p));
Sbeta_bp = aux'*aux/B;

% GLS: Random effect (Sigma is known)
Pond   = Sigma^(-1/2);
Xtilde = kron(eye(N),Pond)*X;
ytilde = kron(eye(N),Pond)*y;
[t_re,St_re,V_re,bhat_re,S_re,D_re]=residual_boots(ytilde,Xtilde,N,'n_replic',B);
Sbeta_re = inv(Xtilde'*Xtilde);
aux      = bsxfun(@minus,t_re,mean(t_re));
Sbeta_bre = aux'*aux/B;

% Computing Kernels
[fp,Fp,bgrid_p]=mykernel(t_p,'kernel',2);
[fre,Fre,bgrid_re]=mykernel(t_re,'kernel',2);

[X1p,X2p] = meshgrid(bgrid_p(:,1),bgrid_p(:,2));
[X1re,X2re] = meshgrid(bgrid_re(:,1),bgrid_re(:,2));

% Printing results
fprintf('\n============================================================\n');
fprintf(' [III] Residual Bootstrap: Panel data, POLS, and GLS\n');
fprintf('============================================================\n');
[bhat_p diag(Sbeta_p).^(1/2) diag(Sbeta_bp).^(1/2) bhat_re diag(Sbeta_re).^(1/2) diag(Sbeta_bre).^(1/2)]
[Sbeta_p(1,2) Sbeta_bp(1,2) Sbeta_re(1,2) Sbeta_bre(1,2)]

% Plots
ww=1;
figure(cc)
h=surf(X1p,X2p,fp'); axis tight;
set(h,'edgealpha',.5,'facelighting','flat','facealpha',.75);
view(-35,30); 
label_x('$$\hat{\beta}_1$$'); label_y('$$\hat{\beta}_2$$');
formataxis(gca); 
% imprime(['Ex4_' num2str(ww)]);
% imprpdf(['Ex4_' num2str(ww)]);
cc=cc+1; ww=ww+1;

figure(cc); hold on;
h1 = plot(t_p(:,1),t_p(:,2));
contour(X1p,X2p,fp',15);
set(h1,'markersize',8,'marker','x','color','blue','linestyle','none');
label_x('$$\hat{\beta}_1$$'); label_y('$$\hat{\beta}_2$$');
xlim([3.5 5.5]); ylim([9 11.5]);
formataxis(gca); 
legend('Data','Isolines'); formatlegend('best'); 
% imprime(['Ex4_' num2str(ww)]);
% imprpdf(['Ex4_' num2str(ww)]);
cc=cc+1; ww=ww+1;

figure(cc)
h=surf(X1re,X2re,fre'); axis tight;
set(h,'edgealpha',.5,'facelighting','flat','facealpha',.75);
view(-35,30); 
label_x('$$\hat{\beta}_1$$'); label_y('$$\hat{\beta}_2$$');
formataxis(gca); 
% imprime(['Ex4_' num2str(ww)]);
% imprpdf(['Ex4_' num2str(ww)]);
cc=cc+1; ww=ww+1;

figure(cc); hold on;
h1 = plot(t_re(:,1),t_re(:,2));
contour(X1re,X2re,fre',15);
set(h1,'markersize',8,'marker','x','color','blue','linestyle','none');
label_x('$$\hat{\beta}_1$$'); label_y('$$\hat{\beta}_2$$');
xlim([3.5 5.5]); ylim([9 11.5]);
formataxis(gca);
legend('Data','Isolines'); formatlegend('best');
% imprime(['Ex4_' num2str(ww)]);
% imprpdf(['Ex4_' num2str(ww)]);
cc=cc+1; ww=ww+1;

%% [V] Pairwise bootstrap (Panel data, continuation)
rng(1115);
% POLS
tpair_p = pairwise_boots(y,X,N,'n_replic',B);
aux      = bsxfun(@minus,tpair_p,mean(tpair_p));
Sbetapair_bp = aux'*aux/B;

% GLS: Random effect (Sigma is known)
tpair_re = pairwise_boots(ytilde,Xtilde,N,'n_replic',B);
aux      = bsxfun(@minus,tpair_re,mean(tpair_re));
Sbetapair_bre = aux'*aux/B;

% Computing Kernels
[fp,~,bgrid_p]=mykernel(tpair_p,'kernel',2);
[fre,~,bgrid_re]=mykernel(tpair_re,'kernel',2);

[X1p,X2p] = meshgrid(bgrid_p(:,1),bgrid_p(:,2));
[X1re,X2re] = meshgrid(bgrid_re(:,1),bgrid_re(:,2));

% Printing results
fprintf('\n============================================================\n');
fprintf(' [IV] Pairwise Bootstrap: Panel data, POLS, and GLS\n');
fprintf('============================================================\n');
[diag(Sbeta_bp).^(1/2)  diag(Sbetapair_bp).^(1/2) diag(Sbeta_bre).^(1/2) diag(Sbetapair_bre).^(1/2)]
[Sbeta_bp(1,2) Sbetapair_bp(1,2) Sbeta_bre(1,2) Sbetapair_bre(1,2)]

% Plots
ww=1;
figure(cc)
h=surf(X1p,X2p,fp'); axis tight;
set(h,'edgealpha',.5,'facelighting','flat','facealpha',.75);
view(-35,30); 
label_x('$$\hat{\beta}_1$$'); label_y('$$\hat{\beta}_2$$');
formataxis(gca); 
% imprime(['Ex5_' num2str(ww)]);
% imprpdf(['Ex5_' num2str(ww)]);
cc=cc+1; ww=ww+1;

figure(cc); hold on;
h1 = plot(tpair_p(:,1),tpair_p(:,2));
contour(X1p,X2p,fp',15);
set(h1,'markersize',8,'marker','x','color','blue','linestyle','none');
label_x('$$\hat{\beta}_1$$'); label_y('$$\hat{\beta}_2$$');
xlim([3.5 5.5]); ylim([9 11.5]);
formataxis(gca); 
legend('Data','Isolines'); formatlegend('best'); 
% imprime(['Ex5_' num2str(ww)]);
% imprpdf(['Ex5_' num2str(ww)]);
cc=cc+1; ww=ww+1;

figure(cc)
h=surf(X1re,X2re,fre'); axis tight;
set(h,'edgealpha',.5,'facelighting','flat','facealpha',.75);
view(-35,30); 
label_x('$$\hat{\beta}_1$$'); label_y('$$\hat{\beta}_2$$');
formataxis(gca); 
% imprime(['Ex5_' num2str(ww)]);
% imprpdf(['Ex5_' num2str(ww)]);
cc=cc+1; ww=ww+1;

figure(cc); hold on;
h1 = plot(tpair_re(:,1),tpair_re(:,2));
contour(X1re,X2re,fre',15);
set(h1,'markersize',8,'marker','x','color','blue','linestyle','none');
label_x('$$\hat{\beta}_1$$'); label_y('$$\hat{\beta}_2$$');
xlim([3.5 5.5]); ylim([9 11.5]);
formataxis(gca);
legend('Data','Isolines'); formatlegend('best');
% imprime(['Ex5_' num2str(ww)]);
% imprpdf(['Ex5_' num2str(ww)]);
cc=cc+1; ww=ww+1;

%% [VI] Wild bootstrap (heteroscedasticity)
rng(1115);   % for reproducibility
% x1 : variables observed at t=1
% x2 : variables observed at t=2

% observable sample
v    = ([x1(:,1) x2(:,1)]).*mvnrnd([0; 0],Sigma,N);
y  = X*btrue+vec(v');
Omg  = [2 1;1 2].*(diag(lambda)*Sigma*diag(lambda));
Pond   = Omg^(-1/2);
Xtilde = kron(eye(N),Pond)*X;
ytilde = kron(eye(N),Pond)*y;

% POLS
[tp_r,~,~,bhat_p,~, Rp]= residual_boots(y, ...
                            X,N,'n_replic',B);
aux   = bsxfun(@minus,tp_r,mean(tp_r));
Sbp_r = aux'*aux/B;

tp_p  = pairwise_boots(y,X,N,'n_replic',B);
aux   = bsxfun(@minus,tp_p,mean(tp_p));
Sbp_p = aux'*aux/B;

tp_w = wild_boots(y,X,N,'n_replic',B,'algo',1);
aux   = bsxfun(@minus,tp_w,mean(tp_w));
Sbp_w = aux'*aux/B;

% GLS estimator
[tg_r,~,~,bhat_g] = residual_boots(ytilde, ...
                            Xtilde,N,'n_replic',B);
Rgp   = reshape(y-X*bhat_g,2,N)';
aux   = bsxfun(@minus,tg_r,mean(tg_r));
Sbg_r = aux'*aux/B;

tg_p  = pairwise_boots(ytilde,Xtilde,N,'n_replic',B);
aux   = bsxfun(@minus,tg_p,mean(tg_p));
Sbg_p = aux'*aux/B;

tg_w  = wild_boots(ytilde,Xtilde,N,'n_replic',B,'algo',1);
aux   = bsxfun(@minus,tg_w,mean(tg_w));
Sbg_w = aux'*aux/B;

% asymptotic estimators
Xr    = bsxfun(@times,x1,Rp(:,1)) + bsxfun(@times,x2,Rp(:,2));
Sa_p  = (X'*X)\((Xr'*Xr)/(X'*X));

Xa = permute(reshape(Xtilde',2,2,N),[3 1 2]);
Rgp2 = reshape(kron(eye(N),Pond)*vec(Rgp'),2,N)';

Xgr   = bsxfun(@times,Xa(:,:,1),Rgp2(:,1)) + bsxfun(@times,Xa(:,:,2),Rgp2(:,2));
Sa_g  = (Xtilde'*Xtilde)\((Xgr'*Xgr)/(Xtilde'*Xtilde));

% Printing results
fprintf('============================================================\n');
fprintf('[V] Wild Bootstrap: Panel data, POLS, and GLS\n');
fprintf('============================================================\n');
[bhat_p diag(Sa_p).^.5 diag(Sbp_r).^.5 diag(Sbp_p).^.5 diag(Sbp_w).^.5; ...
   bhat_g  diag(Sa_g).^.5   diag(Sbg_r).^.5 diag(Sbg_p).^.5 diag(Sbg_w).^.5]
[Sa_p(1,2) Sbp_r(1,2) Sbp_p(1,2) Sbp_w(1,2); Sa_g(1,2) Sbg_r(1,2) Sbg_p(1,2) Sbg_w(1,2)]

% Plots: scatter
ww=1;
figure(cc);
h1 = plot(x1(:,1), abs(Rp(:,1)));
set(h1,'markersize',8,'marker','x','color','blue','linestyle','none');
label_x('var1: t=1'); 
label_y('mod(residuals): t=1'); 
formataxis(gca); 
% imprime(['Ex6_' num2str(ww)]); 
% imprpdf(['Ex6_' num2str(ww)]);
cc=cc+1; ww=ww+1;

figure(cc);
h1 = plot(x2(:,1), abs(Rp(:,2)));
set(h1,'markersize',8,'marker','x','color','blue','linestyle','none');
label_x('var1: t=2'); 
label_y('mod(residuals): t=2');
formataxis(gca); 
% imprime(['Ex6_' num2str(ww)]); 
% imprpdf(['Ex6_' num2str(ww)]);
cc=cc+1; ww=ww+1;

figure(cc);
h1 = plot(x1(:,1), abs(Rgp(:,1)));
set(h1,'markersize',8,'marker','x','color','blue','linestyle','none');
label_x('var1: t=1'); 
label_y('mod(residuals): t=1'); 
formataxis(gca); 
% imprime(['Ex6_' num2str(ww)]); 
% imprpdf(['Ex6_' num2str(ww)]);
cc=cc+1; ww=ww+1;

figure(cc);
h1 = plot(x2(:,1), abs(Rgp(:,2)));
set(h1,'markersize',8,'marker','x','color','blue','linestyle','none');
label_x('var1: t=2'); 
label_y('mod(residuals): t=2');
formataxis(gca); 
% imprime(['Ex6_' num2str(ww)]); 
% imprpdf(['Ex6_' num2str(ww)]);
cc=cc+1; ww=ww+1;

%% [VII] VAR Bootstrap
rng(1115);          % seed for reproducibility    
load dataVAR.mat;   % DD contains data used in VAR
p = 2;              % lag order    
B = 25;            % bootstrap replications [500]
K = 50;             % impulse responses horizon

% Bootstrapping: wild bootstrap
[bhat,Omg,Z,X,b_est,S]=VAR_boots(DD,'order',p,'n_replic',B,'burnin',20,'wild');
[T,m]=size(DD);

% Initialization (media) [Long-run restriction]
BET      = reshape(b_est,m,[]);
LR       = inv(eye(m)-BET(:,2:m+1)-BET(:,m+2:2*m+1));
RA(m^2,m*(m+1)/2) = -LR(4,4)/LR(4,5); 
gamma0 = 0.8*ones(m*(m+1)/2,1); 
gamma0(m*(m+1)/2) = -gamma0(m*(m+1)/2); 
gamma0 = mlSVAR(RA,gamma0,S,T,'tau',0.3, ...
            'errotol_score',1e-7,'errortol_grad',1e-7);
F  = [BET(:,2:end); ...
        [kron(eye(p-1,p-1),eye(m,m)), zeros(m*(p-1),m)]];
A  = reshape(RA*gamma0,m,m);
IR = VARimpulse(F,A,K);

% Bootstrap for IRF's [SVAR]
IRb = [];   
for b=1:B
    beta  = bhat(:,b);
    BET   = reshape(beta,m,[]);
    Sb    = Omg(:,:,b);        
    F = [BET(:,2:end);...
        [kron(eye(p-1,p-1),eye(m,m)), zeros(m*(p-1),m)]];
    eigen = eig(F);          
    if ~any(abs(eigen)>1)   % is stable?
        LR    =  inv(eye(m)-BET(:,2:m+1)-BET(:,m+2:2*m+1));
        % Long-run restriction
        RA(m^2,m*(m+1)/2) = -LR(4,4)/LR(4,5); 
        gammaAb = mlSVAR(RA,gamma0,Sb,T,'tau',.3, ...
                'errotol_score',1e-4,'errortol_grad',1e-6,'quiet');    
        Ab     = reshape(RA*gammaAb,m,m);        
        IRb(:,:,:,b) = VARimpulse(F,Ab,K);
    end
end

% [VII] Plotting IRFs based on structural identification
Kmax = 46;
x  = (1:Kmax)';
ll= 0.16;
uu= 0.84;

sh     = [5 4];             % i / Ds
resp   = [5 4 2 1];         % i/Ds/pi/y
for ss=1:numel(sh)
    for vv=1:numel(resp)
        posit = 'NorthEast';
        tempb=squeeze(IRb(resp(vv),:,sh(ss),:));
        % normalizing
        aux = ones(K+1,1)*squeeze(IRb(sh(ss),1,sh(ss),:))';
        tempb=tempb./aux;
        irf = IR(resp(vv),:,sh(ss))'/IR(sh(ss),1,sh(ss));
        if resp(vv)~=4;
            irfb   = quantile(tempb,.5,2);
            irfb_l = quantile(tempb,ll,2);
            irfb_u = quantile(tempb,uu,2);
        else
            irf=cumsum(irf);
            irfb   = quantile(cumsum(tempb),.5,2);
            irfb_l = quantile(cumsum(tempb),ll,2);
            irfb_u = quantile(cumsum(tempb),uu,2);
        end
        irf=irf(x);
        irfb=irfb(x); irfb_l=irfb_l(x); irfb_u=irfb_u(x);
        miny=min(irfb_l);
        maxy=max(irfb_u); 
        
        figure(cc); hold on;
        dif = irfb_u-irfb_l;
        a = area([irfb_l dif],-15);
        h=plot(x, [irfb irf]);
        plot([1 K],[0 0],':k');
        set(a(1), 'EdgeColor', 'none', 'FaceColor', 'none');
        set(a(2), 'EdgeColor', 'none', 'FaceColor', [0.65 0.65 1]);
        set(h(1), 'Linewidth', 2.5, 'Color', [0.1 0.1 0.65]); 
         set(h(2), 'Linewidth', 2.5, 'Color', [0.65 0.1 0.1],'linestyle','--'); 
        % axis
        ylim([min(miny,0)*1.2 maxy*1.25]);
        set(gca,'xtick',1:5:Kmax,'Xticklabel',0:5:Kmax-1);
        label_x('Months'); label_y('Percent');
        formataxis(gca);
        xlim([0 Kmax+1]);
        if resp(vv)==5
            % legend only for one subplot
            legend([h;a(2)],'IRF median [b]','IRF central',['Percentiles [' num2str(ll) ' - ' num2str(uu) ']']);        
            formatlegend(posit);
        end
%         imprime(['IRF' num2str(resp(vv)) '_' num2str(sh(ss))]);
%         imprpdf(['IRF' num2str(resp(vv)) '_' num2str(sh(ss))]);
        cc=cc+1;
    end
end

%% [VIII] Moving-Block bootstrap
rng(1115);      % for reproducibility
B     = 500;
T     = 200;
p     = 1;
mu    = 2;    
sigma = 2;
rho   = 0.4;
y     = nan(T,1);
y(1)  = 0;
for t=2:T
    y(t)=mu*(1-rho)+rho*y(t-1)+ sigma*randn;
end

lambda = T^1/3;

% MB bootstrap and naive bootstrap
statMB  = MB_boots(y,lambda,@mean,'n_replic',B)';
statNai = naive_boots(y,@mean,'n_replic',B);

% Printing results
rr=[sigma/(1-rho) std([sqrt(T)*(statNai-mean(y)) sqrt(T)*(statMB-mean(y))])];
fprintf('\n==================================\n');
fprintf('[VII] Stationary MB bootstrap \n');
fprintf('\n==================================\n');
fprintf('True  Naive boots  M-B boots\n');
fprintf('==================================\n');
fprintf('%2.4f   %2.4f    %2.4f \n',rr);
fprintf('==================================\n');


toc




