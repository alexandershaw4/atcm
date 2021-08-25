function [v,sPf] = optimise_smoothing(Pf,y,k)

f = @(k) sum( denan((spm_vec(y) - spm_vec(atcm.fun.tsmovavg(Pf','t',k))).^2) );

v = fminsearch(f,k);

sPf = f(v);

