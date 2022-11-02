function [Q,GL] = AGenQn(x,n)
% Convert vector into Gaussian process, i.e. a smoothed symmetric matrix form.
%
% Q = AGenQn(x,n)
%
% Implements function 
%       f = @(x,n) diag(ones(N-abs(x),1),x)/n;
% for the sequence
%       M = f(0,1) + f(1,2) + f(-1,2) + f(2,4) + f(-2,4) + f(3,8) + f(-3,8) ...
%
% AS2022

if nargin < 2 || isempty(n)
    n = 3;
end

N = length(x);

f = @(x,n) diag(ones(N-abs(x),1),x)/n;

% a la, M = f(0,1) + f(1,2) + f(-1,2) + f(2,4) + f(-2,4) + f(3,8) + f(-3,8);

M = f(0,1);

for i  = 1:n
    s  = cumprod(repmat(2,[1 i]));
    M  = M + f(i,s(end)) + f(-i,s(end));
end

% linear model
%Q = smooth2(denan(M.\diag(x)),4);
Q = smooth2(M*diag(x),4);

if nargout == 2
    %Q  = cdist(x,x) .* Q;
    A  = Q .* ~eye(length(Q));
    N  = size(A,1);
    GL = speye(N,N) + (A - spdiags(sum(A,2),0,N,N))/4;
end
