function [t,Y,F,Xgrid,f]=nonparam_boots(X,T,varargin)
% Obtains bootstrap sample for any statistic of interes using
% non-parametric smoothed bootstrap technique.
% ==============
%  Syntax:
% ==============
%       [t,Y,F,Xgrid]=naive_boots(X,T)
%       [t,Y,F,Xgrid]=naive_boots(X,T,...)
% 
% ==============
%   Inputs
% ==============
%  X:   is the observed variables matrix (size: nobs x m)
%  T:   is a handle function for statistic computation.
%
%   **********
%     Options
%   **********
%   [1] 'n_replic': scalar that sets the number of bootstrap samples (B).
%                   [default: B=100]
%   [2] 'n_grid'  : scalar that sets the number of grid in the computation
%                   empirical distribution function. [default: n_grids=50]
%   [3]  'kernel'  : scalar that sets the kernel function used in the kernel 
%                   distribution estimator. [default: Uniform pdf].
%                   See 'mykernel.m' documentation.
%
% ==============
%  Outputs
% ==============
%  t     : bootstrap sample for the statistic of interes.  
%  Y     : bootstrap sample of observed data (based on X).
%  F     : Kernel distribution function (KCE).
%  Xgrid : Grid used in KCE computation.
%  f     : Kernel density estimator (KDE).
%
% ========================================================================
%   By Alex Carrasco Martinez (alex.carmar93@gmail.com), december 2017
% ========================================================================

% [I] Options and settings
B       = 100;
n_grids = 100;
kernel  = 1;
[N,k]   = size(X);

for w=1:numel(varargin)
    if strcmp(varargin{w},'n_replic'), B    = varargin{w+1}; end
    if strcmp(varargin{w},'n_grid'), n_grids = varargin{w+1}; end  
    if strcmp(varargin{w},'kernel'), kernel  = varargin{w+1}; end  
end

% [II] Kernel distribution estimation
[f,F,Xgrid] = mykernel(X,'n_grid',n_grids,'kernel',kernel);

ind_str='ind1';
for m=2:k
    ind_str = [ind_str, [', ind' num2str(m) ] ];
end

% [III] Bootstrapping
Fb = repmat(F,[ones(1,k),N]);
gridsize = size(F);
if k==1
    gridsize = n_grids;
end

Y = nan(N,k,B);
t = [];

for b=1:B  
    rb       = rand(N,1);    
    rb       = reshape(rb,[ones(1,k),N]);
    rb       = repmat(rb,[gridsize,1]); 
    
    whoismin = reshape(Fb-rb,[n_grids^k,N]);       
    whoismin(whoismin>0) = nan; 
    [~,idx]  = max(whoismin);
    eval(['[' ind_str ']=ind2sub(gridsize,transpose(idx));']);    
    eval(['newind = sub2ind([n_grids,k],[' ind_str '], repmat(1:k,N,1));']);
    Y(:,:,b)   = Xgrid(newind);                                      
    t(b,:)     = T(Xgrid(newind));
end

end