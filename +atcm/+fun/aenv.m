function e = aenv(y,n)
% Smoothing for power spectra using envelopes, without shifting peaks
%
% AS2021

S = [min(y) max(y)];

% compute envelope
for i = 1:n
    y = rescale( real(y + sqrt( (y.^2) + (hilbert(gradient(y)).^2)) ),S(1),S(2));
end

e = y;
