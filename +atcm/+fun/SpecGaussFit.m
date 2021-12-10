function model = SpecGaussFit(w,y,n)
% Feature selection function for fitting power spectra - a Gaussian mixture
% model of dim n (n<=8) which returns the model and the model cefficients,
% e.g.
%
% model = atcm.fun.SpecGaussFit(w,y,n)
% model = { m(w) log(coefficients) }
%
%

warning off;
name = sprintf('Gauss%d',n);
m    = fit(w,y,name);
v    = coeffvalues(m);
lv   = log(v);
warning on;

lv(isinf(lv))=0;

model = {m(w) lv};