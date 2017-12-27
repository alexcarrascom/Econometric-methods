function [beta,Omg,Z,X,bhat,S,Vsample]=VAR_boots(D,varargin)
% Obtains parameter estimator resamples using parametric bootstrap for VAR
% models.
%
%    z_t = C_0 + C_1*z_{t-1} + ... + C_p*z_{t-p} + v_t 
%
% where z_t is a (m x 1) random vector, v_t is (m x 1), and C_i is the (m x m) 
% parameter vector for i=1, ..., p.
% ==============
%  Syntax:
% ==============
%       [beta,Omg,Z,X,bhat,S,Vsample]=VAR_boots(D)
%       [beta,Omg,Z,X,bhat,S,Vsample]=VAR_boots(D,...)
% 
% ==============
%   Inputs
% ==============
%  D : Observed sammple of endogenous variables.
% 
%   **********
%    Options
%   **********
%   [1] 'n_replic'   : scalar that sets the number of bootstrap samples (B).
%                      [default: B=100]
%   [2] 'order'      : scalar that sets the number of lag to be considered.
%                      [default: p=1]
%   [3] 'noconstant' : indicates to not include constant in estimation.
%                      [default: const=1]
%   [4] 'wild'       : indicates to use wild bootstrap.
%                      [default: uses residual bootstrap]
%   [5] 'burnin'     : sets how many observations will not be considered in 
%                      parameter estimation.
%                      [default: burn=0]
% ==============
%  Outputs
% ==============
%  beta: (B x k) matrix with vector parameter estimations using bootstrap sample.  
%  Omg : (m x m x B) object with estimation of residual variance matrix.
%  Z   : Stack data for z_t
%  X   : Stack data for z_{t-1},...,z_{t-p}.
%  bhat: OLS estimation using the given sample.
%  S   : residual variance estimation using the given sample.
%  Vsample  : residuals estimation using OLS.
%
% ========================================================================
%   By Alex Carrasco Martinez (alex.carmar93@gmail.com), december 2017
% ========================================================================

% [I] Options and settings
[T,m]   = size(D); 
B       = 100;
p       = 1;
const   = 1;
wild    = 0;
burn    = 0;

for w=1:numel(varargin)
    if strcmp(varargin{w},'n_replic'),   B=varargin{w+1}; end
    if strcmp(varargin{w},'order'),      p=varargin{w+1}; end
    if strcmp(varargin{w},'noconstant'), const=0; end
    if strcmp(varargin{w},'wild'),       wild=1; end
    if strcmp(varargin{w},'burnin'),     burn=varargin{w+1}; end
end

% [II] Step 1: Constructing data
Z = vec(D(p+1:T,:)');
X = [];

for tt=p:T-1
    Zt = vec(D(tt-(0:p-1),:)');
    X = [X,Zt];
end

if const
    Xd=ones(1,T-p);
    X=[Xd;X];    
end

X    = kron(X',eye(m));

% [III] Step 2: Computing estimation
bhat    = (X'*X)\(X'*Z);
v       = Z-X*bhat; 
Vsample = reshape(v,m,T-p)';
S       = Vsample'*Vsample/T;  % ML estimator

% [II] Bootstrapping
TT   = T+burn; 
beta = nan(numel(bhat),B);
Omg  = nan(m,m,B);

for b=1:B
    rs   = randi(T-p,1,TT);
    u    = ones(1,TT);
    if wild
        u = binornd(1,0.5,[1,TT]);
        u = u - (u==0);
    end
    aux  = X(1:m,:);
    Zb   = []; 
    Xb   = aux;      
    for tt=1:TT-1
        y   = aux*bhat+ u(tt)*Vsample(rs(tt),:)';
        aux = [aux(:,1:const*m) kron(y',eye(m)) aux(:,const*m+1:end-m^2)];
        Xb  = [Xb;aux];
        Zb  = [Zb;y];
    end
    y    = aux*bhat+Vsample(rs(TT),:)';
    Zb   = [Zb;y];
    Xb   = Xb(burn*m+1:end,:);
    Zb   = Zb(burn*m+1:end);
    beta(:,b) = (Xb'*Xb)\(Xb'*Zb);
    v    = Zb-Xb*beta(:,b); 
    v    = reshape(v,m,[])';
    Omg(:,:,b)  = v'*v/T;  % ML estimator           
end

end