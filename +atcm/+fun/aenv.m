function e = aenv(y,n,S)
% Smoothing for power spectra using envelopes, without shifting peaks
%
% AS2021

y = real(y);
if nargin < 3 || isempty(S)
    S = [min(y) max(y)];
end

% compute envelope
for i = 1:n
    y = rescale(real( real(y + sqrt( (y.^2) + (hilbert(gradient(y)).^2)) )),S(1),S(2));
end

e = y;
