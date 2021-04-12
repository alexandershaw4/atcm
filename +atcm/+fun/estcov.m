function Cov = estcov(X,MM)

N = length(X);
Y=zeros(N-MM+1,MM);
for m=1:MM
    Y(:,m) = X((1:N-MM+1)+m-1);
end;
Cov=Y'*Y / (N-MM+1);

