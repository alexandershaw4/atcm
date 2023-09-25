function [Ppf,w] = fftgc(pc,dt,w,k)
% compute the fourier transform under Gaussian constraints; for signal pc
% sampled at 1/dt and at frequencies in w;
%
%  [Ppf,w] = atcm.fun.fftgc(pc,dt,w)
%
% by convolving a Gaussian Process matrix with the DFT matrix (basis set)
%
% AS2023

if nargin < 4 || isempty(k)
    k= 18;
end

N = length(pc);
F = dftmtx(N);
G = VtoGauss(ones(size(F)),k,[],0);
F = real(F).*G + sqrt(-1)*(imag(F).*G);
f = (1/dt) * (0:(N/2))/N;

data   = pc*F;
data   = (data/N);
L2     = floor(N/2);
data   = data(1:L2+1);

fwin = find(f<(w(end)+1));
f    = f(fwin);
Ppf  = abs(data(fwin));
Ppf  = interp1(f,Ppf,w,'linear','extrap');

end