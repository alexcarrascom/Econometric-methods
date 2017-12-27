function [Chat,Shat,F,Theta,C0,C,IC]=VARest(D,p,varargin)
% Estimates a reduced form VAR model by OLS, i.e., z_t = C_0 + C_1 z_{t-1}
% + ... C_p z_{t-p}+v_t.
% Syntax:
%      [Chat,Shat,F,Theta,C0,C,IC]=VARest(D,p,varargin)
%
%  [1] D:  (T*m) matrix where m is the number of variables. Note the this
%           matrix is not arranged to obtain OLS estimator directly, this
%           function is designed to make it for us.
%  [2] p:  scalar number that indicates the number of lags to be considered
%          in the estimation.
% Options:
%  [1] 'noconst': Estimate a model without constant.
%  [2] 'determ' : We are able to specfy any deterministic variable to be 
%                 considered in the estimation.
%
% Output:
%   Chat : estimated parameters, [C_1,C_2,..,C_p]. 
%   Shat : covariance matrix for reduced form errors.
%   F    : Matrix of the companion form. (to obtain IRFs)
%   Theta: Long run multipliers
%   C0   : Estimated constant
%   C    : estimated parameters in different arrange (separately).
%   IC   : information criterias 
%
% **********************************************************
%   By Alex Carrasco, november 2017 
% **********************************************************

%% [I] Set-up
cte=1;
determ=0;
T=size(D,1);

for ii=1:numel(varargin)
    if strcmp(varargin{ii},'noconst'), cte=0; end
    if strcmp(varargin{ii},'determ'),  cte=0;
                                       determ=1; 
                                       Xd=varargin{ii+1}; 
    end
end

m=size(D,2);
Z = D(p+1:T,:)';
X = [];

for tt=p:T-1
    Zt = vec(D(tt-(0:p-1),:)');
    X = [X,Zt];
end

% Deterministic components
if cte
    Xd=ones(1,T-p);
    X=[Xd;X];
elseif determ
    X=[Xd;X];
end

%% [II] Estimation
Chat = Z*X'/(X*X');
Shat = (Z-Chat*X)*(Z-Chat*X)'/T;

if cte || determ
    C0   = Chat(:,size(Xd,1));
else
    C0   = [];
end

Theta  = eye(m,m);
C      = nan(m,m,p);
n_d    = size(Xd,1);

for jj=1:p
    m0 =(jj-1)*m+1+n_d;
    C(:,:,jj)=Chat(:,m0:m0+m-1);
    Theta=Theta-C(:,:,jj);
end

Theta=inv(Theta);
F = [Chat(:,1+n_d:end);[kron(eye(p-1,p-1),eye(m,m)), zeros(m*(p-1),m)]]; % companion form

%% [III] information criterias (based on, Lutkhepol (2005))
IC(1)= ((T+p*m+1)/(T-p*m-1))^m*det(Shat);             % FPE (final predictor error)
IC(2)= log(det(Shat)) + p*m^2*2/T;                    % AIC
IC(3)= log(det(Shat)) + p*m^2*log(T)/T;               % BIC
IC(4)= log(det(Shat)) + p*m^2*2*log(log(T))/T;        % HQ

end