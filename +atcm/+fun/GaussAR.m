function [V,K,Ppf] = GaussAR(x,w,dt,N)
% A Gaussian component-based auto-regressive spectral model.
%
%     Pf = atcm.fun.GaussAR(x,w,dt)
%
% Input: x:  timeseries
%        w:  frequencies of interest
%        dt: 1./sample frequency
%
% Process:
% 1 - Uses the Yule-Walker AR method of spectral estimation
% 2 - Conversion to a Gaussian Process, i.e. 2d symmetric smooth matrix
% 3 - PCA of matrix with max limit to 8 components
% 4 - Ensure eigenvectors conform to semi-Gaussian
% 5 - Project back on data
%
% AS

Ppf = pyulear(x(:),8,w,1/dt);

% Make sure it's a mv Gauss so that the error is smooth
Ppm = AGenQn(Ppf);
Ppm = Ppm' ;

% principal components
[V,D] = eig(Ppm);
[~,order] = sort(diag(D),'descend');
D = diag(D);
D = D(order);
V = V(:,order);

% Max components in spectrum
MaxG = 8;

g  = find( (D*pi) > 1 );
ng = length(g);
ng = min(ng,MaxG);
g  = g(1:ng);
V  = V(:,g);

% Ensure eigenvectors are semi-Gaussian
for i = 1:ng
    V(:,i) = AGenQ(V(:,i))*V(:,i);
end

K   = Ppm;

if ng > 1
    Ppf = spm_vec( max((V'*K)) );
else
    Ppf = V'*K;
end

if nargin > 3
    V = V(:,1:N);
end

return;

% Test Gauss-AR
dt = 1/300;
t  = (0:dt:1-dt)';
x  = 4 * sin(2*pi*10*t) + 4 * sin(2*pi*50*t);
w  = 3:80;
[V,~,Pf] = atcm.fun.GaussAR(x(:),w,dt);
figure,plot(w,V);

