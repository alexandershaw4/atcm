function [CC,C,Q] = clustervec(Y,N,q)

if nargin < 3 || isempty(q)
    q = 20;
end
if nargin < 2 || isempty(N)
    N = 5;
end

Q   = atcm.fun.AGenQn(Y,q);
idx = kmeans(Q,N);

for i = 1:N
    C(i,:)  = Y*0;
    fi      = find(idx==i);
    C(i,fi) = Y(fi);
end

for i = 1:N
    CC(i,:) = C(i,:)*atcm.fun.QtoGauss(C(i,:),4);
end

% normalise to Y
CC = CC./sum(CC(:));
CC = CC * sum(Y);

end