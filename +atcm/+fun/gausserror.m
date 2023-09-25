function e = gausserror(Y,y)

dgY = VtoGauss(real(Y));
dgy = VtoGauss(real(y));

Dg  = dgY - dgy;
e   = norm(Dg*Dg') ;