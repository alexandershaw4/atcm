function [W,C] = aregress(X,Y,method,p1,p2)
% Linear model of X on Y;
%
%   [W,C] = atcm.fun.aregress(X,Y,method,p1,p2)
%
%  X(n,m) and Y(n,1)
%  method can be 'ML','MAP' or 'Bayesian'
%  p1 is regulariser on I; p2 for Bayesian method
%
%  Used in the context of gradient descent, where X = partial gradients,
%  the MAP solution is equivalent to a Levenberg-Marquart step.
%
% adapted from github
% AS2024

PHI = X;
D = size(X,2);

if nargin < 4; p1 = 1; end
if nargin < 5; p2 = 1/8; end

if strcmp(method,'ML')==1
    W=(PHI'*PHI)\(PHI'*Y); 
elseif strcmp(method,'MAP')==1
    I=p1.*eye(D);
    W=(I+PHI'*PHI)\(PHI'*Y); 

    % H = I+PHI'*PHI;
    %[E,D] = eig(H);
    %D     = diag(D);
    %H = E*diag(D.*(D < 0))*E';
    % W = H\(PHI'*Y); 

else
    I=eye(D);
    lambda=(p1.^2)./(p2.^2);
    W=(lambda*I+PHI'*PHI)\(PHI'*Y);
    C=inv(p1.^(-2)*(PHI'*PHI)+p2.^(-2).*I);
end

end