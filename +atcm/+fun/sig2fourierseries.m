function [sigout,f,F,b] = sig2fourierseries(sig,k,fs)
% fit signal as k-component Fourier series;
%
%      [sigout,hz,fftmat] = atcm.fun.sig2fourierseries(sig,k)
%
% AS2023

sig = sig(:);

F = atcm.fun.afftmtx(length(sig)*2,k);
F = F(1:length(sig),:);
b = F\sig;

sigout = (F*b);

if nargin < 3 || isempty(fs)
    fs = 1;
end

L  = length(sig);
f  = fs * (0:(L/2))/L;
f  = f(1:k);
