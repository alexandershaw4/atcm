% script to prep and fit a single channel ERP using the TC model

% preparation of the SPM file ERP
%--------------------------------------------------------------------------
Dfile = '/Users/Alex/Dropbox/KET_SPM_VS_Visual/SPM8_VS_KET_POST_02101351_Sens.mat';
D = spm_eeg_load(Dfile);

% condition indices
trials = find( strcmp('Undefined',D.conditions) );

% time indices
toi(1) = atcm.fun.findthenearest(0,  D.time);
toi(2) = atcm.fun.findthenearest(1.5,D.time);

% channel ind
chi = 1;

% trial average data
data = detrend( squeeze(mean( D(chi,toi(1):toi(2),trials) ,3)) );
pst  = D.time(toi(1):toi(2));

% set up a single region model
%--------------------------------------------------------------------------
% function handles
M.IS = @atcm.integrate_erp;
M.f  = @atcm.tcm_nomh;

% time, freqs etc.
M.Hz  = 4:80;

% initial states
M.x = zeros(1,8,7);
M.x(:,:,1) = -70;

DCM.M = M;

% neural priors
DCM = atcm.parameters(DCM,1);

% simulus / input bump
R(1) = 0.69; % input bump (stim) delay: exp(R(1))*60
R(2) = .5;   % input bump (stim) size:  exp(R(2))*8
DCM.M.pE.R = [R 0];

% Trial data
U.X = [0];

DCM.M.sim.pst = pst;
DCM.M.sim.dt  = 1./D.fsample;

% check simulator
[y] = feval(DCM.M.IS,DCM.M.pE,DCM.M,U);

% place real data in xY
xY.y  = {data};
xY.Hz = pst;

% turn off most parameters
pC   = DCM.M.pC;
pC   = spm_unvec( spm_vec(pC)*0, pC);
pC.R = pC.R + 1/8;
pC.L = pC.L + 1/8;
pC.a = pC.a + 1/8;
pC.H = pC.H + eye(8)/16;
pC.T = pC.T + 1/16;
DCM.M.pC = pC;

[Qp,Cp,Eh,F] = atcm.optim.spm_nlsi_GN_as(DCM.M,U,xY);