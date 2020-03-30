function s0 = integrator_pst(P,M,U)
% a wrapper on the spectral RK45 integration scheme defined in integrator
% that returns the time-domain outputs only

[y,w,s,~,~,pst] = feval(M.IS,P,M,U);

s  = s{1};
[ns,np,nk,nt] = size(s);

%s0 = squeeze(s{1}(1,:,1,:));

s0 = reshape(s,[ns*np*nk,nt]);