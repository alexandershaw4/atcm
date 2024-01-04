function dy1 = match_energy_signals(y0,y1,I)
% match energy in signal y1 to that of y0 using svd;
%
% dy1 = match_energy_signals(y0,y1,[I])
%
% optionally retaining I/100% of original signal energy.
%
% AS2023

G = @(x) atcm.fun.VtoGauss(x);

YM = G(spm_vec(y0));
[u0,s0,v0] = svd(YM);

if nargin < 3 || isempty(I)
    I = atcm.fun.findthenearest(cumsum(diag(s0))./sum(diag(s0)),.9);
end

ym = G(spm_vec(y1));
[u,s,v] = svd(ym);

ym  = u(:,1:I)*s(1:I,1:I)*v(:,1:I)';
dy1 = diag(ym);