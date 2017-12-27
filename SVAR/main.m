%% SVAR: An application
% Paper: 'Identifying the interdependence between US monetary policy and
%         the stock market', Bjornland and Leitemo (2009)
% 
% *************************************************************************
%  By Alex Carrasco, november 2017
% *************************************************************************
tic
clear all; clc; close all;
boot = 1;                     % bootstrap again?
B    = 50;                   % numbers of replications [2000]
K    = 100;                   % total periods for IRF

% set up for printing figures 
% path_g = [pwd '\Graphs'];
% imprime = @(x) print( gcf, '-depsc2', [path_g filesep x]);
% imprpdf = @(x) eps2pdf( [path_g filesep x '.eps']);
formataxis  = @(x) set(x, 'Fontname', 'Times', 'FontWeight', 'normal', 'Fontsize', 23, 'Box', 'Off', 'PlotBoxAspectRatio', [1 0.75 1]);
formatlegend  = @(x) set(legend, 'Location', x, 'Orientation', 'Vertical', 'Box', 'Off', 'Fontsize', 22, 'Fontangle', 'normal');

%% [I] loading dataset and data transformation
% [a] loading
load data;
ran = se(1):47:se(end);       % range for plots [do not change!]
m=size(DD,2);

% [b] unit root tests
UR=uroottest(DD,'report','ur_test','names',var_names,'alpha',0.05); % creates a table with test results

% [c] plotting data
cc=1;

for ii=1:m
    figure(cc);
    h=plot(se,DD(:,ii),'-','linewidth',2.5); axis tight;   
    set(h,'color',[0 0 0.5]);
    x=gca;
    set(x,'xtick',ran,'xticklabel',ranlab);
    formataxis(x);
%     imprime([cod_names{ii} '_data']);
%     imprpdf([cod_names{ii} '_data']);
    cc=cc+1;
end

%% [II] Model selection, p: FPE, AIC, BIC, HQ
pmax=8;
[popt,ic]=VAROptLag(DD,pmax);
fprintf('\n Optimal lags: %1.0f %1.0f %1.0f %1.0f\n', popt);

%% [III] Reduced form estimation: VAR(p) 
p = popt(3);     % BIC optimal lag

% [a] Baseline model
[Chat,Shat,F,Theta]=VARest(DD,p);

% [b] Alternative model
DD(:,[4 5])=DD(:,[5 4]);
[Chata,Shata,Fa]=VARest(DD,p);

DD(:,[4 5])=DD(:,[5 4]);

P = chol(Shat);                  % Cholesky decomposition (baseline)
Pa = chol(Shata);                % Cholesky decomposition (alternative)           

%% [IV] A-model, identification and estimation
% RA is the restriction matrix
A = nan(m,m);
% [a] Zero restrictions
A(1,2:5) = 0; A(2,3:5)=0; A(3,4:5)=0;

free_param = isnan(A);
nfree      = sum(vec(free_param));

ind0 = reshape(1:m^2,m,m);
ind0 = ind0(free_param);

RA = zeros(m^2,nfree);
[I,J] = ind2sub([m,m],ind0);
for k=1:nfree
    RA((J(k)-1)*m+I(k),k)=1;
end

% [b] Long-run restriction
RA(m^2,m*(m+1)/2)=-Theta(4,4)/Theta(4,5);
RA(:,m*(m+1)/2+1)=[];

% [c] Estimation: 
%  Based on Amisano and Gianini (1997)/ Marquardt algorithm (Baseline ordenation)
gamma0=.8*ones(nfree-1,1); gamma0(nfree-1)=-gamma0(nfree-1); % Sign restriction
[gammaA,logL,error,S]=mlSVAR(RA,gamma0,Shat,T,'tau',0.3,'errotol_score',1e-7,'errortol_grad',1e-7);
A=reshape(RA*gammaA,m,m);

%% [V] IRF [mean]
% [a] Structural identifcation
IR=VARimpulse(F,A,K); % mean responses in structural identifcation

% [b] Models identificated with recursive restrictions
Kmax = 46;        % toplot
x =(1:Kmax)';
sh = [5 4];       % shocks: i/Ds \\ Ds/i
resp = [5 4];     % shocks: i/Ds \\ Ds/i
IRchol1 =VARimpulse(F,P',K);    % IRF using baseline order
IRchol2 =VARimpulse(Fa,Pa',K);  % IRF using alternative order

for ss=1:numel(sh)
    for vv=1:numel(resp)
        figure(cc);
        irf1 = IRchol1(resp(vv),:,sh(ss))'/IRchol1(sh(ss),1,sh(ss));  % normalizing responses
        irf2 = IRchol2(resp(resp(vv)~=resp),:,sh(sh(ss)~=sh))'/IRchol2(sh(sh(ss)~=sh),1,sh(sh(ss)~=sh));
        if resp(vv)==4
            irf1=cumsum(irf1);
            irf2=cumsum(irf2);
        end
        irf1=irf1(x); irf2=irf2(x); 
        
        h=plot([irf1,irf2]); hold on;
        plot([1 Kmax],[0 0],':k');
        set(h(1),'color',[0 0 .5],'linewidth',3);
        set(h(2),'color',[0 0 .5],'linestyle','--','linewidth',2);
        
        % axis
        axs=gca; miny=min(min([irf1,irf2])); maxy=max(max([irf1,irf2])); 
        xlim([1 Kmax]); 
        ylim([min(miny,0)*1.2 maxy*1.2]);
        set(axs,'xtick',6:5:Kmax,'Xticklabel',5:5:Kmax-1);
        xlabel('Months','fontname','Times','FontWeight', 'normal', 'Fontsize', 22); 
        ylabel('Percent','fontname','Times','FontWeight', 'normal', 'Fontsize', 22);
        formataxis(axs);
        if resp(vv)==5
            lg=legend('Baseline','Alterntive');        
            formatlegend('NorthEast');
        end
            
%         imprime(['IRFChol_' cod_names{resp(vv)} '_' cod_names{sh(ss)}]);
%         imprpdf(['IRFChol_' cod_names{resp(vv)} '_' cod_names{sh(ss)}]);
        cc=cc+1;
    end
end

%% [VI] Standard errors [Bootstraping]
if boot                             
    [bhat,Omg]=VAR_boots(DD,'order',p,'n_replic',B,'burnin',20,'wild');
    [T,m]=size(DD);
    
    fprintf('\n Estimating bootstrap sample ... '); kk=1;    
    IRb = [];   
    for b=1:B
        if b==floor(kk*B/10)
            fprintf('%1.0f ',kk);
            kk=kk+1;
        end
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
    fprintf(', done!\n');

    save SVAR_irf IR IRb DD var_names B
else
    load SVAR_irf
end

%% [VII] Plotting IRFs based on structural identification
Kmax = 46;
x  = (1:Kmax)';
ll=0.16;
uu=0.84;

sh = [5 4];             % i / Ds
resp   = [5 4 2 1];     % i/Ds/pi/y
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
            irfb = quantile(cumsum(tempb),.5,2);
            irfb_l = quantile(cumsum(tempb),ll,2);
            irfb_u = quantile(cumsum(tempb),uu,2);
        end
        irf = irf(x);
        irfb=irfb(x); irfb_l=irfb_l(x); irfb_u=irfb_u(x);
        miny=min(irfb_l);
        maxy=max(irfb_u); 
        
        figure(cc); hold on;
        dif = irfb_u-irfb_l;
        a = area([irfb_l dif],-15);
        h=plot(x, irfb);
        plot([1 K],[0 0],':k');
        set(a(1), 'EdgeColor', 'none', 'FaceColor', 'none');
        set(a(2), 'EdgeColor', 'none', 'FaceColor', [0.65 0.65 1]);
        set(h, 'Linewidth', 2.5, 'Color', [0 0 0.5]); 
        % axis
        axs=gca;     
        ylim([min(miny,0)*1.2 maxy*1.2]);
        set(axs,'xtick',1:5:Kmax,'Xticklabel',0:5:Kmax-1);
        xlabel('Months','fontname','Times','FontWeight', 'normal', 'Fontsize', 22); 
        ylabel('Percent','fontname','Times','FontWeight', 'normal', 'Fontsize', 22);
        formataxis(axs);
        xlim([0 Kmax+1]);
        if resp(vv)==5
            % legend only for one subplot
            legend([h;a(2)],'IRF',['Percentiles [' num2str(ll) ' - ' num2str(uu) ']']);        
            formatlegend(posit);
        end
%         imprime(['IRF' cod_names{resp(vv)} '_' cod_names{sh(ss)}]);
%         imprpdf(['IRF' cod_names{resp(vv)} '_' cod_names{sh(ss)}]);
        cc=cc+1;
    end
end
toc
%% Extra [Model selection table]
fprintf('\n');
for ii=1:size(ic,1)
fprintf('$%2.0f$ & $%2.4f$ & $%2.4f$ & $%2.4f$ & $%2.4f$ \\\\ \n',ii,ic(ii,:))
end