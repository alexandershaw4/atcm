function [y,s,u,Q] = gaulinsvdfit(x,n)

if nargin < 2
    n = 1;
end

Q = atcm.fun.VtoGauss(x);
m = diag(Q) - x(:);
for i = 1:length(Q)
    Q(i,:) = Q(i,:) - m(i);
end

[u,s,v] = svd(Q - x);
s = diag(s);
y = spm_vec(s(1:n)'*u(:,1:n)'*Q);

s = s(1:n)';
u = u(:,1:n)';
Q = Q;