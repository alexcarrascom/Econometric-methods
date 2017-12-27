function [popt,IC]=VAROptLag(DD,pmax,varargin)
% Selects optimal lag in a VAR model based on information criterias.
% popt: is a vector with the optimal lag for FPE, AIC, BIC , HQ (in that order)
% ************************************
%  By Alex Carrasco, november 2017
% ************************************

IC=nan(pmax,4);
for p=1:pmax
    [~,~,~,~,~,~,ica]=VARest(DD,p);
    IC(p,:)=ica;
end
[~,popt]=min(IC);

end