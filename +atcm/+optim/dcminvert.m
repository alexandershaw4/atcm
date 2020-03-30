function DCM = dcminvert(DCM)
% use the default DCM Bayeisan EM+ inversion routine, & save


switch char(DCM.M.IS)
    case 'atcm.integrate_kern'
        
        % use the (fast) kernel inversion
        DCM.M.IS     = 'spm_csd_mtf';
        [Qp,Cp,Eh,F] = atcm.optim.spm_nlsi_GN_as(DCM.M,DCM.xU,DCM.xY);
        % then switch to full integrator for full predictions
        DCM.M.IS     = @atcm.integrate_kern;
        Ns           = length(DCM.A{1});
        
    otherwise

        % Variational Laplace: model inversion
        %==========================================================================
        [Qp,Cp,Eh,F] = atcm.optim.spm_nlsi_GN_as(DCM.M,DCM.xU,DCM.xY);
        Ns           = length(DCM.A{1});
end

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
 
% unpack priors for Bayeian inference
%--------------------------------------------------------------------------
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

% Call the model with full outputs, using posterior parameters
%--------------------------------------------------------------------------
[Hs Hz dtf] = feval(DCM.M.IS,qp,M,DCM.xU);
[~,~,s,ModSig,InputSig,pst,layers,noisemod,firing,DelayMat] = feval(DCM.M.IS,qp,DCM.M,DCM.xU);
dt = (pst(2)-pst(1))/1000;
[ccf]       = spm_csd2ccf(Hs,DCM.M.Hz,dt);
[coh fsd]   = spm_csd2coh(Hs,DCM.M.Hz);
DCM.dtf     = dtf;
DCM.ccf     = ccf;
DCM.coh     = coh;
DCM.fsd     = fsd;

% Store the model staimates
%--------------------------------------------------------------------------
DCM.pst                 = pst;
DCM.Hz                  = Hz;
DCM.AllModelStates      = s;
DCM.ModelWeightedOutput = ModSig;
DCM.InputSig            = InputSig;
DCM.LayerSpectra        = layers;
DCM.NoiseComponents     = noisemod;
DCM.FiringP             = firing;   % note: convert to Hz using atcm.fun.f2sr
DCM.Delays              = DelayMat;
 
% Compute functional amplitude & phase coupling between all cells (all regions)
%--------------------------------------------------------------------------
DCM.FunctionalConnectivty = atcm.fun.computefc(s);
[DCM.phases,DCM.PhsCor]   = atcm.fun.computephase(s);

% store parameter estimates in DCM
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