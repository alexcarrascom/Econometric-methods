function [gammaA,logL,error,S0]=mlSVAR(RA,gamma0,SigmaV,T,varargin)
% mlSVAR: Computes quasi maximum likelihood estimator for an structural VAR 
% with 'model A' residual structure whenever identification restrictions
% can be expressed as: RA*gammaA=vec(A).
% 
% [gamma,logL,error,S0] = mlSVAR(R,gamma0,sigmaV,T,options)
% Lavenberg-Marqueardt Newton - Raphson method is used as numeric
% algorithm.
% ***************************************************
%   by Alex Carrasco, october 2017
% ***************************************************

m=size(SigmaV,1);
errortol1 = 1e-7;
errortol2 = 1e-3;
maxiter  = 1e5;  
iter=1;
k=1; % step
tau=0;
quiet = 0;
logL = [];

for i = 1:length(varargin)
     if strcmpi(varargin{i}, 'maxiter')          ; maxiter = varargin{i+1}; end
     if strcmpi(varargin{i}, 'errortol_grad')    ; errortol1 = varargin{i+1}; end
     if strcmpi(varargin{i}, 'errortol_score')   ; errortol2 = varargin{i+1}; end
     if strcmpi(varargin{i}, 'kstep')            ; k = varargin{i+1}; end
     if strcmpi(varargin{i}, 'tau')              ; tau = varargin{i+1}; end
     if strcmpi(varargin{i}, 'quiet')            ; quiet = 1; end
end

logF = @(iA) -m*T*log(2*pi)/2 + T*log(det(iA))-0.5*T*trace(iA'*iA*SigmaV);
A    = reshape(RA*gamma0,m,m);
iA   = inv(A);
logL0 = logF(iA);
K    = commutation(m,m);
iSigmaV = inv(SigmaV);
aux2 = .5*T*kron(iSigmaV,iSigmaV);
S    = -RA'*kron(iA,iA')*(T*vec(A')- T*kron(SigmaV,eye(m))*vec(iA));
aux1 = (eye(m^2)+K)*kron(A,eye(m));
I    = RA'*aux1'*aux2*aux1*RA;
mu   = tau*max(diag(I));
v    = 2;
nfree = numel(gamma0);

if ~quiet
    fprintf('\n F.I.M.L estimation for a SVAR model: A - model\n');
    fprintf('===============================================================\n');
    d={'Nº Iter' 'log Likelihood' 'Score' 'Changes'};
    fprintf('%s %2s %s %3s %s %4s %s\n',d{1},'',d{2},'',d{3},'',d{4});
end

stop=0;
while ~stop && (iter <= maxiter) %|| (norm(S)>errortol)    
    if rem(iter,500)==0 && ~quiet;
        fprintf('%5.0f %4s %3.4f %3s %2.4f %6s %2.4f\n',iter,'',logL0,'' ,norm(S),'',norm(Dg));
    end
    rho=-1;
    while (rho<0) && ~stop 
        N = I+mu*eye(nfree,nfree);
        Dg = k*(N\S);
        gamma1 = gamma0 + Dg;
        A = reshape(RA*gamma1,m,m);
        iA = inv(A);
        logL(iter) = logF(iA);
        DL = logL(iter)-logL0;
        
        rho = DL/(Dg'*(mu*Dg+S));
        if (norm(Dg)<= errortol1)
            stop=1;
        elseif rho>0
            stop = (norm(S)<= errortol2) || (norm(logL0)<= errortol2); 
            S    = -RA'*kron(iA,iA')*(T*vec(A')- T*kron(SigmaV,eye(m))*vec(iA));
            aux1 = (eye(m^2)+K)*kron(A,eye(m));
            I    = RA'*aux1'*aux2*aux1*RA;
            gamma0 = gamma1;
            logL0 = logL(iter);
            mu=mu*max([1/3, 1-(2*rho-1)^2]); v=2;
        else
            mu = mu*v; v=2*v;
        end
    end
    iter=iter+1;
end 
if ~quiet
    fprintf('%5.0f %4s %3.4f %3s %2.4f %6s %2.4f\n',iter,'',logL0,'',norm(S), '',abs(norm(Dg)));
end

gammaA=gamma0;
error = norm([Dg;DL]);
S0 = S;

end