function [DCM] = convert_to_timefreq(DCM)

hz = [10 80];
bl = [-1000 0];

% augment params
DCM.M.pE.U_delay = 0;
DCM.M.pC.U_delay = 1/8;

DCM.M.pE.C_delay = 0;
DCM.M.pC.C_delay = 1/8;


P = DCM.M.pE;
M = DCM.M;
U = DCM.xU;

% compute data timefreq
[TFR, f, t, TFR_trials] = atcm.tf.dcm_data_tfm(DCM,hz,bl);


% make vg mask
W = atcm.tf.make_visual_gamma_mask(t, f);

data.TFR = TFR;
data.W = W;

tfopts.fband = hz;
tfopts.baseline = bl;

M.fs = 1./DCM.xY.dt;
M.T = t;



% compute model tf
[y_vec, cache] = atcm.tf.tcm_tfr_forward_vectorised(P, M, U, tfopts, data);


% make function we can optimise
f = @(p) atcm.tf.tcm_tfr_forward_vectorised(spm_unvec(p,P), M, U, tfopts, data);

p0 = spm_vec(DCM.M.pE);
S0 = diag(spm_vec(DCM.M.pC));

[m, V, D, logL, iter, sigma2, allm,g_elbo] = ...
    fitVariationalLaplaceThermo(TFR(:)*0, f, p0, S0, 32, 1e-6,1);

Dinv  = spdiags(1./diag(D), 0, size(D,1), size(D,1));
Mid   = eye(size(V,2)) + V'*(Dinv*V);        % k×k
CP    = Dinv - Dinv*V*(Mid \ (V'*Dinv));

DCM.Ep = spm_unvec(m,DCM.M.pE);
DCM.Cp = CP;