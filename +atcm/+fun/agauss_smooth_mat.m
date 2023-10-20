function [P,ym,uu] = agauss_smooth_mat(x,n)
% this function performs sliding-window smoothing with an n-width Gaussian 
% kernel on the input vector x (of length k), stacking each window into a matrix of
% size kxk (a sort of Gaussian-process-like matrix), then perofrms PCA on 
% this matrix and returns the PCA components explaining 99% of cumulative sum
% of eigenvalues.
%
% [P,matrix] = agauss_smooth_mat(x,n)
%
% AS2023

x   = x(:);
w   = (1:length(x))';
fun = @(Wid,f) 2 * exp( -(w-f).^2 / (2*(2*Wid)^2) );

for i = 1:length(x)
    ym(i,:) = ( x.*fun(n,i) );
end

% decompose matrix using svd of covariance
[u,s,v] = svd(cov(ym'));

i = atcm.fun.findthenearest(cumsum(diag(s))./sum(diag(s)),.99);
P = u(:,1:i)'*ym';
P = abs(P);
uu = u(:,1:i);

end