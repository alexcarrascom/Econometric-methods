function IR=VARimpulse(F,A,K)
% F: companion form
% A: identificated matrix for structural schoks
% K: Impulse response horizon
% IR: structural impulse responses
% *****************************************
%   By Alex Carrasco, november 2017
% *****************************************

m = size(A,1);
IR=[];
for k=0:K
	aux = F^k;
    IR(:,k+1,:) = aux(1:m,1:m)*A;
end

end