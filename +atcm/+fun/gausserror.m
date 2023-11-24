function e = gausserror(Y,y)

Y = spm_vec(Y);
y = spm_vec(y);


dgY  = VtoGauss(real((Y)));
dgy  = VtoGauss(real((y)));
Dg   = (dgY - dgy).^2;

e = ((dgY - dgy).^2)/Y(:)';

e = norm(e,'fro');

%e    = norm(Dg*iS*Dg','fro');

% peaks?
p0  = atcm.fun.indicesofpeaks(real(Y));
p1  = atcm.fun.indicesofpeaks(real(y));
dp  = cdist(p0(:),p1(:));
if isvector(dp)
    dp = abs(diag(dp));
end

dp = denan(dp);

peake = trace(diag(diag(dp)));

peake = denan(peake);
peake = abs(peake);
peake = max(peake,1/2);

e   = abs(e) * abs(peake);


end

% 
% 
% 
% dgY = VtoGauss(real(Y));
% dgy = VtoGauss(real(y));
% 
% Dg  = dgY - dgy;
% e    = norm(Dg*Dg','fro');
% 
% % peaks?
% p0  = atcm.fun.indicesofpeaks(real(Y));
% p1  = atcm.fun.indicesofpeaks(real(y));
% dp = cdist(p0(:),p1(:));
% if isvector(dp)
%     dp = diag(dp);
% end
% 
% e   = e * trace(diag(diag(dp)));