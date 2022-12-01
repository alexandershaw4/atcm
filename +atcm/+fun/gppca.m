function y = gppca(Ppf)
% a sort of gaussian process decomposition of a signal or matrix
%
% y = gppca(Ppf)
%
% AS


% comute covariance and smooth with Gauss kernel
if all(size(Ppf)>1)
    Ppm = cov(Ppf');
    Ppm = Ppm.*atcm.fun.GaussEye(size(Ppm,2))*Ppm';
else
    MM = atcm.fun.estcov(Ppf,length(Ppf));
    Ppm = sqrt(MM*atcm.fun.GaussEye(length(Ppf))*MM');
end

% eigenvectors are spectral modes
[V,D] = eig(Ppm);
[~,order] = sort(diag(D),'descend');
D = diag(D);
D = D(order);
V = V(:,order);

Ppf = diag(V*V');%*X;
y   = AGenQ(Ppf)*Ppf;