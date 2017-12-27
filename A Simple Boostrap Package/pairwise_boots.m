function [beta,Yb,Xb]=pairwise_boots(y,X,nobs,varargin)
% Obtains parameter estimator resamples using pairwise bootstrap method 
% in a linear model:
%
%    y_i = X_i*beta + v_i
%
% where y_i is a (m x 1) random vector, X_i is (m x k), v_i is (m x 1), 
% and beta is the (k x 1) parameter vector.
% ==============
%  Syntax:
% ==============
%       [beta,Yb,Xb]=residual_boots(y,X,nobs)
%       [beta,Yb,Xb]=residual_boots(y,X,nobs,...)
% 
% ==============
%   Inputs
% ==============
%  y:   is a stack vector of y_i for i=1 to nobs (size: m.nobs x 1 )
%  X:   is a stack vector of X_i for i=1 to nobs (size: m.nobs x k )
%  nobs: is the number of observations in sample
%   **********
%     Options
%   **********
%       'n_replic': scalar that sets the number of bootstrap samples (B).
%                   [default: B=100]
% ==============
%  Outputs
% ==============
%  beta : (B x k) matrix with vector parameter estimations using bootstrap sample.  
%  Yb   : (nobs.m x B) object with dependent variable bootstrap sample.
%  Xb   : (nobs.m x k x B) object with covariates bootstrap sample.
%
% ========================================================================
%   By Alex Carrasco Martinez (alex.carmar93@gmail.com), december 2017
% ========================================================================

% [I] Options
B       = 100;
m       = numel(y)/nobs;

for w=1:numel(varargin)
    if strcmp(varargin{w},'n_replic'), B=varargin{w+1}; end
end

% [II] Bootstrapping
[~,k] = size(X);
Ydat  = reshape(y,m,nobs)';
Xdat  = reshape(X',k,m,nobs);
Yb    = nan(nobs*m,B);
Xb    = nan(nobs*m,k,B);
beta  = nan(B,k);

for b=1:B  
    rs           = unidrnd(nobs,1,nobs);
    Yb(:,b)      = vec(Ydat(rs,:)');
    auxX         = permute(Xdat(:,:,rs),[2 3 1]);
    Xb(:,:,b)    = squeeze(reshape(auxX,[nobs*m,k,1]));
    beta(b,:)    = (Xb(:,:,b)'*Xb(:,:,b))\(Xb(:,:,b)'*Yb(:,b));
end

end