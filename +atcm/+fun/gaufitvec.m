function y = gaufitvec(w,Ppf)

X0 = atcm.fun.gauregmat(w);

Ppf = real(Ppf);

warning off;

b = X0\Ppf;

warning on;

b(b<0)=0;

m = [min(Ppf) max(Ppf)];

Ppf = X0*b;

y = rescale(Ppf,m(1),m(2));