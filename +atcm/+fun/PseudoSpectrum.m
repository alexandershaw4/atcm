function X = PseudoSpectrum(w,f0)
% Generate a pseudo power spectrum at freqs w with peaks at f0
%
% AS

n = length(f0);
X = atcm.fun.makef(w,f0,f0./w(end),ones(1,n)*4).*(w.^-2);
X = X.*rescale(hamming(length(X))');