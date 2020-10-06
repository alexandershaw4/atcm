function [EP,F,CP,History,M] = default_fit(DCM)


M = AODCM(DCM);

M.opts.maxit        = 4;
M.opts.hyperparams  = 1;
M.opts.BTLineSearch = 0;

M.optimise();

% posteriors->priors for a re-run using smaller steps
M.update_parameters(M.Ep);
M.opts.step_method = 3;

M.optimise();

% And a final re-run again using bigger steps
M.update_parameters(M.Ep);
M.opts.step_method = 1;

M.optimise();


% Extract the things we need
EP = M.Ep;
F  = M.F;
CP = M.CP;
History = M.history;

EP = spm_unvec( spm_vec(EP), DCM.M.pE);