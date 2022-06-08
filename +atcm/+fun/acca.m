function [A,B] = acca(X,Y)

% Center the variables
X = X - mean(X,1);
Y = Y - mean(Y,1);

[Q1,T11,perm1] = qr(X,0);
[Q2,T22,perm2] = qr(Y,0);

[n,p1] = size(X);
rankX = rank(X);
rankY = rank(Y);

d = min(rankX,rankY);
[L,D,M] = svd(Q1' * Q2,0);
A = T11 \ L(:,1:d) * sqrt(n-1);
B = T22 \ M(:,1:d) * sqrt(n-1);
r = min(max(diag(D(:,1:d))', 0), 1); % remove roundoff errs

% Put coefficients back to their full size and their correct order
A(perm1,:) = A;
B(perm2,:) = B;
