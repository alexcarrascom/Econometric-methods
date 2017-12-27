function [f,F,Xgrid]=mykernel(X,varargin)
% Computes kernel density estimator (KDE) and kernel cumulative estimator (KCE). 
% ==============
%  Syntax:
% ==============
%       [f,F,Xgrid]=mykernel(X)
%       [f,F,Xgrid]=mykernel(X,...)
% 
% ==============
%   Inputs
% ==============
%  X:   is the observed variables matrix (size: nobs x m)
%   **********
%     Options
%   **********
%       'kernel'  : scalar that sets the kernel function used in the kernel 
%                   distribution estimator:
%                   [1] kernel=1, Uniform density
%                   [2] kernel=2, Gaussian density
%                   [2] kernel=3, Epanechnikov kernel function
%                   [4] kernel=4, Triangular kernel function
%                   [default: kernel=1].
%       'n_grid'  : scalar that sets the number of grid in the computation
%                   empirical distribution function. [default: n_grids=50]
%       'bandwith': sets the bandwith size. [Default: use Scott's rule of 
%                   thumb, h=(4/((k+2)*N))^(1/(k+4))*std(X)]
%     'eachseries': Compute KDE for each series in database.  
% ==============
%  Outputs
% ==============
%  f     : Kernel density estimator (KDE).
%  F     : Kernel distribution function (KCE).
%  Xgrid : Grid used in KCE computation.
%
% ========================================================================
%   By Alex Carrasco Martinez (alex.carmar93@gmail.com), december 2017
% ========================================================================

% [I] Options
n_grids  = 100;
kernel = 1;
[N,k]  = size(X);
h      = (4/((k+2)*N))^(1/(k+4))*std(X); % Scott's rule of thumb
minX   = min(X)-h; maxX=max(X)+h;
flageach = 0;

for w=1:numel(varargin)
    if strcmp(varargin{w},'kernel'),     kernel=varargin{w+1}; end
    if strcmp(varargin{w},'n_grid'),     n_grids=varargin{w+1}; end
    if strcmp(varargin{w},'bandwidth'),  h=varargin{w+1}; end
    if strcmp(varargin{w},'eachseries'), flageach=1; end
end


% [II] Setting up
Xgrid  = ones(n_grids,1)*minX + diag(0:(n_grids-1))*ones(n_grids,1)*(maxX-minX)/(n_grids-1); 

dom='X1';  gg  = 'Xgrid(:,1)';
for m=2:k
    dom  = [dom, [' X' num2str(m) ] ];
    gg   = [gg, [', Xgrid(:,' num2str(m) ')'] ];
end
eval( ['[' dom ']=ndgrid(' gg ');'] );

switch kernel
   case 1 %'Uniform'
        K  = @(psi) (abs(psi)<=1)*1/2;
        IK = @(psi) (abs(psi)<=1).*(psi+1)/2 + (psi>1);
    case 2 %'Gaussian'
        K  = @normpdf;
        IK = @normcdf;
    case 3 %'Epanechnikov'
        K  = @(psi) (abs(psi)<=1).*(3/4*(1-psi.^2));
        IK = @(psi) (abs(psi)<=1).*(0.5+3/4*psi -1/4*psi.^3) + (psi>1); 
    case 4 %'Triangular'
        K  = @(psi) (abs(psi)<=1).*(1-abs(psi));
        IK = @(psi) (psi>=-1 && psi<0).*((1-abs(psi))*(psi+1)/2) ...
                    +(psi>=0 && psi<1).*(1-(1-abs(psi))*(-psi+1)/2) ...
                    + (psi>=1);
end

gridsize=size(X1);
f = zeros(gridsize);
F = zeros(gridsize);

% [III] Kernel estimator
if ~flageach
    for i=1:N
        aux1=1;
        aux2=1;
        for j=1:k
            eval(['aux3 = (X(i,j)-X' num2str(j) ')/h(j);']);
            aux1 = aux1.*K(aux3);     
            aux2 = aux2.*IK(-aux3);
        end
        f=f+aux1;
        F=F+aux2;
    end
    
    f=f/(N*prod(h));
    F=F/N;
    
else
    for i=1:k
        f(:,i)=ind_kernel(X(:,i),Xgrid(:,i),K,h(i),N,n_grids);
    end
    F=[];
end


end

function y=ind_kernel(x,Xgrid,k,h,N,n_grids)

y=(N*h)^(-1)*sum(k((Xgrid*ones(1,N)-(x*ones(1,n_grids))')/h),2);

end