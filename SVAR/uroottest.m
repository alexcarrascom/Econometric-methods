function UR=uroottest(DD,varargin)
% Executes ADF and PP unit root test with different lar order for univariate
% time series and selects the best model based on BIC.
% Syntax:
%      UR=uroottest(DD,varargin)
%  [1] DD: (T*m) matrix where m is the number of variables.
% Options:
%  [1] 'lag_max': maximum number of lags for selecting the model.
%  [2] 'alpha'  : Size of the test.
%  [3] 'test'   : Cell object specifying one of both tests (ADF/PP).
%  [4] 'H1'     : Cell object specifying the model consistent with the
%                 alternative hypothesis. 'ard' is a stationary AR(p) with
%                 constant and 'ts' is a trend stationary AR(p).
%  [5] 'report' : Make a .tex file with the main results for ADF and PP
%                 test. OBS: you must give a name for the file.
% Output:
%   UR: Double which their rows contains results for the different test. 
% **********************************************************************
%   By Alex Carrasco, November 2017
%
%   This function needs some functions of the Econometric Toolbox
% **********************************************************************

%% [I] Default
report = 0;
m = size(DD,2);
pmax = 8;
alpha = 0.1;
test   = {'adf','pp'};
test_h1 = {'ard','ts'};
names      = cellstr([repmat('var ',m,1),num2str((1:m)')]);
%% [II] Options
for ii=1:numel(varargin)
    if strcmp(varargin{ii},'lag_max'), pmax=varargin{ii+1}; end
    if strcmp(varargin{ii},'alpha'),   alpha=varargin{ii+1}; end
    if strcmp(varargin{ii},'test'),    test=varargin{ii+1}; end
    if strcmp(varargin{ii},'H1'),      test_h1=varargin{ii+1}; end
    if strcmp(varargin{ii},'names'),   names=varargin{ii+1}; end
    if strcmp(varargin{ii},'report'),  txtname=varargin{ii+1}; report=1; end
end

%% [III] Executing 
warning('off');
select = 1:pmax+1;
UR=nan(m,numel(test)*numel(test_h1)*2);
fprintf('\n Unit root testing ... ');
for hh=1:m
    mod=1;
    for ww=1:numel(test)
        for qq=1:numel(test_h1)
            switch test{ww}
            case 'adf'
                [r,~,~,~,reg]=adftest(DD(:,hh),'model',test_h1{qq},'lags',select-1,'alpha',alpha);
            case 'pp'
                [r,~,~,~,reg]=pptest(DD(:,hh),'model',test_h1{qq},'lags',select-1,'alpha',alpha);
            end
            reject = any(r);
            UR(hh,2*mod-1)=reject;
            candim=select(r==reject);  % select the best model among rejected model or not rejected models
            bestm = candim(end);
            candim=rot90(candim(1:end-1),2);
            for rr=candim
                [~,ind]=min([reg(rr).BIC,reg(bestm).BIC]);
                bestm=bestm+(rr-bestm)*(ind==1);
            end
            UR(hh,2*mod)=bestm-1;      % obtaining the effective number of lags
            mod=mod+1;
        end
    end    
end
fprintf('done!\n');

%% [IV] text and reports
if report
    fid = fopen([txtname '.tex'] ,'wt');

    fprintf(fid, ['\\hspace*{+0.2cm} \\begin{tabular}{l'  repmat('c',1,numel(test)*numel(test_h1)*2) '}\n']);
    fprintf(fid,'\\toprule \n');
    % Header
    fprintf(fid,'\\multirow{3}[3]{*}{\\textsc{\\textbf{Variable}}}');
    for ii=1:numel(test)
        fprintf(fid,'& \\multicolumn{%2.0f}{c}{\\textsc{\\textbf{%s}}}',numel(test_h1)*2,test{ii});
    end
    fprintf(fid,'\\\\ \n');
    % Lines 1
    for ii=1:numel(test)
        fprintf(fid,'\\cmidrule(lr){%1.0f-%1.0f}',(ii-1)*numel(test_h1)*2+2,ii*numel(test_h1)*2+1);
    end
    % Different models
    for ii=1:numel(test)
        for mm=1:numel(test_h1)
            fprintf(fid,'& \\multicolumn{2}{c}{\\textsf{%s}}',test_h1{mm});
        end
    end
    fprintf(fid,'\\\\ \n');
    % Lines 2
    mod=1;
    for ii=1:numel(test)
        for mm=1:numel(test_h1)
            fprintf(fid,'\\cmidrule(lr){%1.0f-%1.0f}',mod+1,mod+2);
            mod=mod+2;
        end
    end
    fprintf(fid,[' ' repmat('& \\textit{Reject \\textsc{H}_0?} & \\textit{Lag} ',1,numel(test)*numel(test_h1))]);
    fprintf(fid,'\\\\ \n');
    fprintf(fid,'\\midrule \n');

    for ii=1:m
        temp=names{ii};
        fprintf(fid,['\\textit{%s} ' repmat(' & %2.0f',1,numel(test)*numel(test_h1)*2) ''],temp,UR(ii,:));
        fprintf(fid,'\\\\ \n');
    end
    fprintf(fid,'\\bottomrule \n');
    fprintf(fid, '\\end{tabular} \\hspace*{+0.2cm} \n');
end


end