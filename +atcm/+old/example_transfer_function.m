
% load a DCM
load old30nov23/aFP_TCM_NewMean_020415_50_control.mat

% change integrator handle to transfer function
DCM.M.IS = @atcm.fun.alex_tf;

% e.g. [y,w] = feval(DCM.M.IS,DCM.M.pE,DCM.M,DCM.xU);

% find fixed point initial states;
x = atcm.fun.alexfixed(DCM.M.pE,DCM.M);

DCM.M.x = spm_unvec(x,DCM.M.x);

