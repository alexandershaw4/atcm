function auto = csd2auto(csd)
% auto = csd2auto(csd)
% assumes csd(nf,ns,ns)
% AS

[nf,ns,ns] = size(csd);
ccsd       = reshape(csd,nf,ns.^2,[]);
auto       = ccsd(:,1:(ns+1):end);