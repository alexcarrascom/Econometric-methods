function [statboot,Dboot]=MB_boots(D,lambda,stat,varargin)
% Executes moving-block bootstrap using D [tobs x m matrix].
% ==============
%  Syntax:
% ==============
%       [statboot,Dboot]=MB_boots(D,lambda,stat)
%       [statboot,Dboot]=MB_boots(D,lambda,stat,...)
% 
% ==============
%   Inputs
% ==============
%  D      : Observed sammple.
%  lambda : mean of the size for each block.
%  stat   : function handle with statistic formula.
%
%   **********
%     Options
%   **********
%   [1] 'n_replic'  : scalar that sets the number of bootstrap samples.
%                    [default: B=100]
%   [2] 'size'      : scalar that sets the size of each resample.
%                    [default: T=tobs].
% ==============
%  Outputs
% ==============
%  statboot: (k x B) matrix with bootstrap sample of the statistic of interest.  
%  Dboot   : (T x m x B) object with resamples.
%
% ========================================================================
%   By Alex Carrasco Martinez (alex.carmar93@gmail.com), december 2017
% ========================================================================

% [I] Options and settings
[tobs,m]  = size(D); 
B         = 100;
T         = tobs; 

for w=1:numel(varargin)
    if strcmp(varargin{w},'n_replic'),   B=varargin{w+1}; end
    if strcmp(varargin{w},'size'),       T=varargin{w+1}; end
end

p=1/lambda;

% [II] Bootstrapping
Dboot    = nan(T,m,B);
statboot = [];

for b=1:B
    count=0;
    X    = [];
    while count<T
        li = geornd(p);
        ti = randi(T);
        try
            X = [X;D(ti:ti+li-1,:)];
            count=count+li; 
        catch
            X = [X;D(ti:end,:)];
            count=count+T-ti+1;
        end        
    end
    X=X(1:T,:);
    statboot(:,b)   = stat(X);
    Dboot(:,:,b)    = X;
end

end