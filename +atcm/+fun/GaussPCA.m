function J = GaussPCA(X,N)

if nargin < 2 || isempty(N)
    N = 1;
end

for i = 1:size(X,2)

    QM = atcm.fun.AGenQn(X(:,i),8);
    [u,s,v] = svd(QM);
    J(i,:) = sum( QM*v(:,1:N), 2);

end

J = abs(J);

end