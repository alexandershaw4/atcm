function [f,J,D] = tcm_delay(x,u,P,M)

M.f = @atcm.tc_hilge2;
[v,J,D] = atcm.tc_hilge2(x,u,P,M);

dt = M.sim.dt;

D     = real(D);
D_dt  = (D*1000)*dt;

b    = pinv(full(J)'.*x(:)).*v(:);
Q    = J.*b; % dxdt = Q*x; Q is linear operator
f    = (Q + (Q.*D_dt) )*v;

