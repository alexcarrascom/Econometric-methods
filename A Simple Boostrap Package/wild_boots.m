function [beta,V,bhat,S,D]=wild_boots(y,X,nobs,varargin)
% Obtains parameter estimator resamples using wild bootstrap method 
% in a linear model:
%
%    y_i = X_i*beta + v_i
%
% where y_i is a (m x 1) random vector, X_i is (m x k), v_i is (m x 1), 
% and beta is the (k x 1) parameter vector.
% ==============
%  Syntax:
% ==============
%    [beta,V,bhat,D]=wild_boots(y,X,nobs)
%    [beta,V,bhat,D]=wild_boots(y,X,nobs,...)
% 
% ==============
%   Inputs
% ==============
%  y:   is a stack vector of y_i for i=1 to nobs (size: m.nobs x 1 )
%  X:   is a stack vector of X_i for i=1 to nobs (size: m.nobs x k )
%  nobs: is the number of observations in sample
%
%   **********
%    Options
%   **********
%   [1] 'n_replic'  : scalar that sets the number of bootstrap samples (B).
%                   [default: B=100].
%   [2] 'size'      : scalar that sets the size of each bootstrap sample.
%                   [default: N=nobs].
%   [3] 'algorithm' : scalar (1, 2, or 3) that sets the algorithm used to
%                     compute auxiliar variable. 1, implements
%                     Rademacher distribution. 2, implements Mammmen's two 
%                     point distribution. 3, uses a Gaussian distibution.
%                    [default: 1].
%   [4] 'HC'        : logic (0 or 1) which specifies residual transformation 
%                     function. 0, does not make any transformation. 1,
%                     uses f(e)=(1-h)^0.5*e. [default: 0]
%                       
% ==============
%  Outputs
% ==============
%  beta : (B x k) matrix with vector parameter estimations using bootstrap sample.  
%  V   : (nobs x m x B) object with residual bootstrap sample (transformed).
%  bhat: OLS estimation using the given sample.
%  S   : residual variance estimation using the given sample.
%  D   : residuals estimation using OLS.
%
% =========================================================================
%   By Alex Carrasco Martinez (alex.carmar93@gmail.com), december 2017
% =========================================================================

% [I] Options
B       = 100;
m       = numel(y)/nobs;
algo    = 1;
hc      = 0;
N       = nobs;

for w=1:numel(varargin)
    if strcmp(varargin{w},'n_replic'), B=varargin{w+1}; end
    if strcmp(varargin{w},'algorithm'),algo=varargin{w+1}; end
    if strcmp(varargin{w},'HC'),       hc=varargin{w+1}; end
    if strcmp(varargin{w},'size'),     N=varargin{w+1}; end
end

switch algo
    case 1 %'F2' 
         p=1/2; aux1=1; aux2=-1; 
    case 2 %'F1'
         p=(sqrt(5)+1)/(2*sqrt(5)); aux1=-(sqrt(5)-1)/2; aux2=(sqrt(5)+1)/2;   
end

% [II] Step 1
[~,k]  = size(X);          % k is the number of parameters
bhat   = (X'*X)\(X'*y);
Px     = X*((X'*X)\X');
v      = y-X*bhat; 
D      = reshape(v,m,nobs)';
S     = D'*D/(nobs-k);

if hc
    for i=1:nobs
        from = (i-1)*m+1;
        to   =   from+m-1;
        D(i,:) = (eye(m)-Px(from:to,from:to))^(-1)*D(i,:)';
    end
end

% [II] Step 2: Bootstrapping
V     = nan(N,m,B);
beta  = nan(B,k);

for b=1:B  
    if algo==3
        u  = randn(N,1);
    else
        u  = aux1*binornd(1,p,[N,1]);
        u  = u + (u==0)*aux2;
    end
    V(:,:,b)  = bsxfun(@times,D,u);
    beta(b,:) = bhat + (X'*X)\(X'*vec(V(:,:,b)'));
end

end