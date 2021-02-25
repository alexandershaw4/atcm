function [y,w] = FastFreqInt(P,M,U)
%the eigensolution frequency-space integrator

w     = M.Hz;                     % FoI (w)
x     = M.x;                      % model (hidden) states
M.x   = spm_dcm_neural_x(P,M);

[f,dfdx,D] = feval(M.f,M.x,M.u,P,M);

dfdx  = D*dfdx;
[v,s] = eig(full(dfdx),'nobalance');
s     = diag(s);

dfdu  = spm_diff(M.f,M.x,M.u,P,M,2);
dfdu  = D*dfdu;
dvdu  = pinv(v)*dfdu;


% use makef to position Xhz oscillations in frequency space 
S = atcm.fun.makef(w,s,dvdu,s*0+4);
 

% for k = 1:length(s)
%     Sk = 1./(1j*2*pi*w - s(k));
%     S  = S + dvdu(k,:)*Sk;
% end

Pf = S(:);

% [Gu,Gs,Gn] = spm_csd_mtf_gu(P,M.Hz);
% 
% Pf = (Pf.*Gs) + Gn;
% 
% for i = 1:length(Pf)
%     Pf(i) = Pf(i)*Gu(i)*Pf(i);
% end

y = {exp(P.L)*Pf};