function [y,b,c,M] = adct(x,target,pct)
% spectral filtering/ smoothing (and fitting) using DCT
%
% examples:
%   xdot = atcm.fun.adct(x);            % dct smoothing
%   y    = atcm.fun.adct(x,target_vec); % dct series linear model of target
%
%   [y,b,c,M] = adct(x,target); % also returns the weights (b) and vectors
%                               (c) such that y = b'*c and c(i,:) = x(i)*M(i,:)
% AS2021

if nargin < 3 || isempty(pct)
    pct = 0.99;
end

x = x(:);
X = dct(x);
[XX,ind] = sort(abs(X),'descend');
i = 1;
while norm(X(ind(1:i)))/norm(X) < pct
  i = i + 1;
end
needed = i;

%needed = atcm.fun.findthenearest( cumsum(abs(X(ind)))./sum(abs(X(ind))), pct );

X(X<(mean(X)-std(X))) = 0;

X(ind(needed+1:end)) = 0;
%X(ind(1))=0;
y = idct(X);

b = 0;
c = [];
%M = dctmtx(length(y));

if nargin > 1 && ~isempty(target)
    % compute weights on x's DCT matrix that produce 'target' 
    
    M = dctmtx(length(y));
    for i = 1:length(y)
        c(i,:) = y(i)*M(i,:);
    end
    
    b = atcm.fun.lsqnonneg(real(c),real(target(:)));
    %b  = inv(c*c')*c*target(:);
    y  = b'*c;
    
end