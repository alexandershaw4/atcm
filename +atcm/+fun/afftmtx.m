function y = afftmtx(N,K)

x = ones(1,N);

for n = 1:N
    for k = 1:K
        y(n,k) = x(n)*exp(-j*2*pi*(k-1)*(n-1)/N);
    end
end