function [f,J,D] = tcm_delay(x,u,P,M)

M.f = @atcm.tc_hilge2;
[v,J,D] = atcm.tc_hilge2(x,u,P,M);

dt = M.sim.dt;

D   = real(D);
Q   = (1 - D*dt);%.*(~~real(J));%inv(eye(npp*nk) - D);
QJ  = Q.*J;

g    = v(:) - x(:);
b    = J'\g;
f    = x(:) + QJ'*b ;