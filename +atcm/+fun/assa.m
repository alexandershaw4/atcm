function [RC,lam] = assa(X,MM,reduce)
% Singular spectrum analysis 
%
%   RC = assa(X,MM);
%
% Optionally returns the n-component projection of the original data, X:
%
%   RC = assa(X,MM,n);
%
% AS

if nargin < 3 || isempty(reduce)
    reduce=0;
end
%MM = 30;                             % delay/embedding
N = length(X);
Y=zeros(N-MM+1,MM);
for m=1:MM
    Y(:,m) = X((1:N-MM+1)+m-1);
end;
Cemb=Y'*Y / (N-MM+1);
C=Cemb;

[RHO,LAMBDA] = eig(C);
LAMBDA = diag(LAMBDA);               % extract the diagonal elements
[LAMBDA,ind]=sort(LAMBDA,'descend'); % sort eigenvalues
RHO = RHO(:,ind);                    % and eigenvectors

PC = Y*RHO;

lam=LAMBDA;

RC=zeros(N,MM);
for m=1:MM
    buf=PC(:,m)*RHO(:,m)'; % invert projection
    buf=buf(end:-1:1,:);
    for n=1:N % anti-diagonal averaging
        RC(n,m)=mean( diag(buf,-(N-MM+1)+n) );
    end
end

% optionally return the n-component projection only
if reduce
    pcx = corr(X,RC).^2;
    [~,I]=sort(pcx,'descend');
    I = I(1:reduce);
    RC = sum(RC(:,I),2);
end