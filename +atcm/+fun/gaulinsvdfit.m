function y = gaulinsvdfit(x)

Q = atcm.fun.VtoGauss(x);
m = diag(Q) - x(:);
for i = 1:length(Q)
    Q(i,:) = Q(i,:) - m(i);
end

[u,s,v] = svd(Q - x);
y = spm_vec(u(:,1)'*Q);