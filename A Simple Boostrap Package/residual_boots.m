function [beta,St,V,bhat,S,D]=residual_boots(y,X,nobs,varargin)
% Obtains parameter estimator resamples using residual bootstrap method 
% in a linear model:
%
%    y_i = X_i*beta + v_i
%
% where y_i is a (m x 1) random vector, X_i is (m x k), v_i is (m x 1), 
% and beta is the (k x 1) parameter vector.
% ==============
%  Syntax:
% ==============
%      [beta,St,V,bhat,S,D]=residual_boots(y,X,nobs)
%      [beta,St,V,bhat,S,D]=residual_boots(y,X,nobs,...)
% 
% ==============
%   Inputs
% ==============
%  y:   is a stack vector of y_i for i=1 to nobs (size: m.nobs x 1 )
%  X:   is a stack vector of X_i for i=1 to nobs (size: m.nobs x k )
%  nobs: is the number of observations in sample.
%   **********
%    Options
%   **********
%   [1] 'n_replic': scalar that sets the number of bootstrap samples (B).
%                   [default: B=100].
%   [2] 'size'   : scalar that sets the size of each bootstrap sample.
%                   [default: N=nobs].
% ==============
%  Outputs
% ==============
%  beta: (B x k) matrix with vector parameter estimations using bootstrap sample.  
%  St  : (m x m x B) object with variance estimations for residuals.
%  V   : (nobs x m x B) object with residual bootstrap sample.
%  bhat: OLS estimation using the given sample.
%  S   : residual variance estimation using the given sample.
%  D   : residuals estimation using OLS.
%
% ========================================================================
%   By Alex Carrasco Martinez (alex.carmar93@gmail.com), december 2017
% ========================================================================

% [I] Options and settings
B       = 100;
m       = numel(y)/nobs;
N       = nobs;

for w=1:numel(varargin)
    if strcmp(varargin{w},'n_replic'), B=varargin{w+1}; end
    if strcmp(varargin{w},'size'),     N=varargin{w+1}; end
end

% [II] Step 1
[~,k] = size(X);          % k is the number of parameters
bhat  = (X'*X)\(X'*y);
v     = y-X*bhat; 
D     = reshape(v,m,nobs)';
S     = D'*D/(nobs-k);

% [III] step 2: Bootstrapping
V    = nan(N,m,B);
beta = nan(B,k);
St   = nan(m,m,B);

for b=1:B      
    rs         = unidrnd(N,1,N);
    V(:,:,b)   = D(rs,:);
    beta(b,:)  = bhat+(X'*X)\(X'*vec(V(:,:,b)'));
    St(:,:,b)  = V(:,:,b)'*V(:,:,b)/(N-k);
end

end