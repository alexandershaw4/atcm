
P = DCM.M.pE;
M = DCM.M;

dt = 1/600;
t  = 0:dt:2-dt;  x0 = DCM.M.x(:);
t0 = 0.0; dur = 2; f = 20;   % 20 Hz for 0.5 s
u  = @(tt) (double(tt >= t0 & tt < t0+dur) .* sin(2*pi*f*(tt - t0)))/128;
opts = struct('method','rosenbrock', 'substeps',8, 'autoAlpha',true, 'finiteCheck',true);
[X, info] = atcm.tcm_integrate_staticJD(DCM.M.f, M, P, x0, t, u, opts);

plot(t,X)