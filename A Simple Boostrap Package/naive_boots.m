function [t,Y,F,Xgrid]=naive_boots(X,T,varargin)
% Obtains bootstrap sample for any statistic of interest using naive 
% bootstrap method. 
% 
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
%   **********
%     Options
%   **********
%       'n_replic': scalar that sets the number of bootstrap samples (B).
%                   [default: B=100]
%       'n_grid'  : scalar that sets the number of grid in the computation
%                   empirical distribution function. [default: n_grids=50]
%       'size'    : scalar that sets the size of each bootstrap sample.
%                   [default: N= nobs]
%
% ==============
%  Outputs
% ==============
%  t     : bootstrap sample for the statistic of interest.  
%  Y     : bootstrap sample of observed data (based on X).
%  F     : empirical distribution function (EDF).
%  Xgrid : Grid used in EDF computation.
%
% ========================================================================
%   By Alex Carrasco Martinez (alex.carmar93@gmail.com), december 2017
% ========================================================================

% [I] Options and settings
B         = 100;
n_grids   = 50;
[nobs,k]  = size(X);
N         = nobs;

for w=1:numel(varargin)
    if strcmp(varargin{w},'n_replic'), B=varargin{w+1}; end
    if strcmp(varargin{w},'n_grid')  , n_grids=varargin{w+1}; end   
    if strcmp(varargin{w},'size')    , N=varargin{w+1}; end   
end

% [II] Empirical distribution function
minX   = min(X); maxX=max(X);
Xgrid  = ones(n_grids,1)*minX + diag(0:(n_grids-1))*ones(n_grids,1)*(maxX-minX)/(n_grids-1); 
F      = nan(n_grids,1);
for i=1:n_grids
    aux = ones(N,1)*Xgrid(i,:);
    F(i)   = sum(X<=aux)/nobs;
end

% [III] Bootstrapping
Y = nan(N,k,B);
t = [];
for b=1:B  
    rs       = unidrnd(N,1,N);
    Y(:,:,b) = X(rs,:);
    t(b,:)     = T(Y(:,:,b));
end

end