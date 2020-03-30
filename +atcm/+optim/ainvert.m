function DCM = ainvert(DCM,method,np)
% non-default inversion using an alternative optimiser & optional parameter
% reduction (clustering)
%
% methods: 'nlsi' - the DCM variational Bayesian routine 
%          'ps'   - particel swarm optimisation
%          'sa'   - simulated annealing
%          'ls'   - laplacian sampling
%          'pattern' - pattern search
%          'hj'   - Hooke-Jeeves pattern
%          'surrogate' - surrogate model optimisation
%          'abc' - artificial bee colony
%
% + lots more - see atcm.optim.clusteroptimse help

if nargin < 3 || isempty(np)
    np = length(find(spm_vec(DCM.M.pC)));
end

% inversion with Alex's parameter-custering routine:
%==========================================================================
[X0,X,cm,F] = atcm.optim.clusteroptimse(DCM.M.pE,DCM,method,np);

Qp = spm_unvec(X0,DCM.M.pE);
Cp = diag(spm_vec(DCM.M.pC));
Ns = length(DCM.A{1});


% Data ID
%--------------------------------------------------------------------------
try
    try
        ID = spm_data_id(feval(DCM.M.FS,DCM.xY.y,DCM.M));
    catch
        ID = spm_data_id(feval(DCM.M.FS,DCM.xY.y));
    end
catch
    ID = spm_data_id(DCM.xY.y);
end
 
pE = DCM.M.pE;
pC = DCM.M.pC;

% Bayesian inference {threshold = prior} NB Prior on A,B and C = exp(0) = 1
%==========================================================================
warning('off','SPM:negativeVariance');
dp  = spm_vec(Qp) - spm_vec(pE);
Pp  = spm_unvec(1 - spm_Ncdf(0,abs(dp),diag(Cp)),Qp);
warning('on', 'SPM:negativeVariance');
 
 
% predictions (csd) and error (sensor space)
%--------------------------------------------------------------------------
Hc  = feval(DCM.M.IS,Qp,DCM.M,DCM.xU);                   % prediction
Ec  = spm_unvec(spm_vec(DCM.xY.y) - spm_vec(Hc),Hc);     % prediction error
 
 
% predictions (source space - cf, a LFP from virtual electrode)
%--------------------------------------------------------------------------
M             = rmfield(DCM.M,'U'); 
M.dipfit.type = 'LFP';

M.U         = 1; 
M.l         = Ns;
qp          = Qp;
qp.L        = ones(1,Ns);             % set virtual electrode gain to unity
qp.b        = qp.b - 32;              % and suppress non-specific and
qp.c        = qp.c - 32;              % specific channel noise

[Hs Hz dtf] = feval(DCM.M.IS,qp,M,DCM.xU);
[~,~,s,ModSig,InputSig,pst,layers,noisemod,firing,DelayMat] = feval(DCM.M.IS,DCM.Ep,DCM.M,DCM.xU);
%[ccf pst]   = spm_csd2ccf(Hs,DCM.M.Hz);
%[coh fsd]   = spm_csd2coh(Hs,DCM.M.Hz);
DCM.dtf     = dtf;
%DCM.ccf     = ccf;
%DCM.coh     = coh;
%DCM.fsd     = fsd;
DCM.pst     = pst;
DCM.Hz      = Hz;
DCM.AllModelStates = s;
DCM.ModelWeightedOutput = ModSig;
DCM.InputSig = InputSig;
DCM.LayerSpectra = layers;
DCM.NoiseComponents = noisemod;
DCM.Firing = firing;
DCM.Delays = DelayMat;

 
% store estimates in DCM
%--------------------------------------------------------------------------
DCM.Ep = Qp;                   % conditional expectation
DCM.Cp = Cp;                   % conditional covariance
DCM.Pp = Pp;                   % conditional probability
DCM.Hc = Hc;                   % conditional responses (y), channel space
DCM.Rc = Ec;                   % conditional residuals (y), channel space
DCM.Hs = Hs;                   % conditional responses (y), source space
DCM.Ce = exp(-Eh);             % ReML error covariance
DCM.F  = F;                    % Laplace log evidence
DCM.ID = ID;                   % data ID
 
% and save
%--------------------------------------------------------------------------
DCM.options.Nmodes = Ns;
 
[fp fn fe] = fileparts(DCM.name); % ensure local saving
if isempty(fp); 
    namename = [fn fe];
else
    namename = [fp '/' fn fe];
end
DCM.name = namename

save(DCM.name, 'DCM', spm_get_defaults('mat.format'));