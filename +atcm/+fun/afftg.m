function [Ppf,w] = afftg(ys,dt,w)

ys  = ys(:);
w   = w(:);
ys  = atcm.fun.bandpassfilter(ys,1/dt,[w(1) w(end)]);

F   = dftmtx(length(ys));
N   = length(F);
f   = (1/dt) * (0:(N/2))/N;
Mt  = pinv(atcm.fun.cdist(f(:),w(:)))';


for i = 1:size(Mt,2)
    %Mt(:,i) = atcm.fun.VtoGauss(Mt(:,i))*Mt(:,i);
    Mt(i,:) = Mt(i,:)*atcm.fun.VtoGauss(Mt(i,:));
end

% Pf(w) = fft(y)*iM, where iM is the inverse distance between
% natural frequency vector and those of interest, subject to
% each column of iM conforming to a Gaussian

Ppf = (ys'*F(:,1:length(f))*Mt)./N;
Ppf = abs(Ppf(:));