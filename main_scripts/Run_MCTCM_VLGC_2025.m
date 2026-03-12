function Run_MCTCM_VLGC_2025(i)
% Run_MCTCM_VLGC_2025
%
% Top-level pipeline for fitting a multi-compartment thalamo-cortical neural
% mass model (MC-TCM) to M/EEG cross-spectral data using a delay-aware
% Laplace-domain transfer function and Variational Laplace optimisation in
% generalised coordinates.
%
% OVERVIEW
% -------------------------------------------------------------------------
% This script implements a Dynamic Causal Modelling (DCM)-style workflow for
% electrophysiological data using an extended thalamo-cortical model that
% augments the standard conductance-based TCM with:
%
%   - separate somatic and dendritic compartments
%   - distinct pre- and postsynaptic contributions
%   - short-term synaptic plasticity (STP) states
%   - compartment-specific routing of excitatory and inhibitory inputs
%
% The aim is to fit M/EEG cross-spectral density (CSD) data with a richer
% circuit model that can capture dendrite-versus-soma dynamics and dynamic
% modulation of synaptic efficacy, while retaining efficient frequency-domain
% inversion through local linearisation around a stable operating point.
%
% In broad terms, the pipeline is:
%
%   1) Load and prepare M/EEG spectral data
%   2) Specify the multi-compartment thalamo-cortical generative model
%   3) Define intrinsic / extrinsic connectivity and exogenous input structure
%   4) Expand and initialise the hidden state space
%   5) Set priors over biophysical, synaptic, compartmental, and STP parameters
%   6) Search for a stable fixed point of the nonlinear system
%   7) Linearise the dynamics around that operating point
%   8) Evaluate a delay-aware transfer function in the Laplace / frequency domain
%   9) Compare predicted and observed spectra
%  10) Optimise model parameters using Variational Laplace in generalised
%      coordinates
%  11) Save fitted posteriors, predictions, and diagnostics into a DCM struct
%
%
% MODEL
% -------------------------------------------------------------------------
% The hidden-state dynamics are defined by the multi-compartment neural mass
% model:
%
%   DCM.M.f = @atcm.tc_twocmp_stp
%
% This model extends the original thalamo-cortical neural mass formalism by
% introducing additional state variables and biophysical structure. The
% hidden states evolve according to:
%
%   dx/dt = f(x,u,P,M)
%
% where:
%
%   x   - hidden neuronal and synaptic states
%   u   - exogenous input
%   P   - model parameters
%   M   - model structure and metadata
%
% Relative to the simpler single-compartment TCM, this version includes:
%
%   - somatic membrane voltage states
%   - dendritic membrane voltage states
%   - presynaptic resource / utilisation variables for STP
%   - compartment-specific synaptic routing
%   - axial soma-dendrite coupling
%
% This makes the model more biophysically expressive and allows the fitted
% spectra to reflect not only synaptic gain and time constants, but also
% dendritic filtering, presynaptic dynamics, and compartment-specific
% inhibition / excitation.
%
%
% MULTI-COMPARTMENT / STP STATE SPACE
% -------------------------------------------------------------------------
% The state space is expanded from the standard 7-state formulation to a
% 10-state formulation per population. In addition to the standard membrane
% and conductance states, the model introduces:
%
%   Vd   - dendritic membrane voltage
%   R    - synaptic resource variable
%   uSTP - dynamic synaptic utilisation / use variable
%
% This creates a richer operating point and allows the model to express
% activity-dependent changes in synaptic transmission. Initial values for
% these additional states are chosen to place the system in a plausible
% baseline regime before fixed-point refinement.
%
%
% DATA PREPARATION
% -------------------------------------------------------------------------
% The script expects a text file listing SPM-format electrophysiology
% datasets. These are read and converted into DCM-compatible cross-spectral
% form.
%
% Preparation steps include:
%
%   - selecting trials / conditions
%   - choosing channels and source labels
%   - restricting the frequency range of interest
%   - estimating / smoothing spectra
%   - assembling the DCM.xY observation structure
%
% The helper function:
%
%   atcm.fun.prepcsd
%
% is used to prepare observed spectra for fitting. The final model inversion
% is performed on CSD data rather than raw time series.
%
%
% TRANSFER FUNCTION FORMULATION
% -------------------------------------------------------------------------
% The forward model is evaluated using a custom delay-aware transfer
% function:
%
%   DCM.M.IS = @atcm.fun.Alex_LaplaceTFwD
%
% Rather than repeatedly simulating the full nonlinear model in the time
% domain during optimisation, this approach:
%
%   - finds a stable operating point
%   - linearises the nonlinear hidden-state equations around that point
%   - constructs the local Jacobian and input mappings
%   - incorporates delays in the frequency domain
%   - computes the model’s spectral response over the frequencies of interest
%
% In local linear form, the dynamics are approximated as:
%
%   dx/dt ≈ A x + B u
%
% and the frequency response is evaluated using a Laplace / resolvent
% formulation of the form:
%
%   H(s) = C * (sI - A_eff(s))^(-1) * B
%
% where A_eff contains delay-dependent phase factors. This yields an
% efficient and numerically tractable way to predict spectra and cross-
% spectra from a nonlinear conductance-based model.
%
% For this multi-compartment model, the transfer-function approximation is
% especially useful because the underlying nonlinear system is larger and
% more expensive to simulate directly.
%
%
% COMPARTMENT-SPECIFIC OBSERVATION MODEL
% -------------------------------------------------------------------------
% The observation model is configured so that different hidden-state
% components can contribute differently to the measured signal. In
% particular, observation weights are constructed using:
%
%   atcm.build_twocmp_observer_J
%
% which allows the measured signal to reflect different contributions from
% somatic and dendritic dynamics across populations. This is important in
% multi-compartment models because the observed field potential is generally
% not a simple readout of a single state, but a weighted combination of
% subcellular and population-level processes.
%
%
% PRIORS AND PARAMETERISATION
% -------------------------------------------------------------------------
% Priors are initialised using:
%
%   DCM = atcm.parameters(DCM,Ns)
%
% and then modified manually to accommodate the extended multi-compartment
% model. The prior structures pE and pC include parameters governing:
%
%   - observation weights (J)
%   - synaptic gains and conductances
%   - delay / smoothing terms
%   - source / channel scaling
%   - axial soma-dendrite coupling (gc)
%   - routing of excitation to dendrites (w_dend)
%   - routing of inhibition to soma (w_soma)
%   - short-term plasticity parameters:
%         prel   baseline release probability
%         tauR   recovery time constant
%         tauU   facilitation / utilisation time constant
%         U0     baseline synaptic use
%
% These priors regularise inference and encourage the model to remain in a
% biologically plausible regime. They also define which aspects of the
% circuit are allowed to vary during inversion.
%
%
% SHORT-TERM PLASTICITY (STP)
% -------------------------------------------------------------------------
% A key feature of this model is the inclusion of short-term synaptic
% plasticity. STP is represented through hidden states and parameters that
% govern moment-to-moment changes in effective synaptic transmission.
%
% In particular, the model includes:
%
%   R    - available synaptic resources
%   uSTP - current synaptic utilisation
%
% with priors on:
%
%   prel - baseline release strength / probability
%   tauR - resource recovery time constant
%   tauU - utilisation / facilitation time constant
%   U0   - baseline use parameter
%
% These mechanisms allow the model to express activity-dependent synaptic
% depression / facilitation and thereby shape spectral output in a way that
% is not possible with a purely static synaptic gain formulation.
%
%
% FIXED POINT / OPERATING POINT SEARCH
% -------------------------------------------------------------------------
% Before spectral fitting, the script searches for a stable fixed point of
% the nonlinear state equations:
%
%   f(x*,0,P,M) = 0
%
% using:
%
%   atcm.fun.find_fixed_point_robust
%
% The initial STP utilisation state is seeded from the prior U0 parameter via
% a logistic transform before fixed-point search, so that the operating point
% is internally consistent with the prior synaptic-use regime.
%
% The fixed point is crucial because the Laplace-domain transfer function is
% based on a local linearisation around this operating point. In effect, the
% fitted model describes the spectral dynamics of the system near its inferred
% baseline state.
%
%
% OPTIMISATION
% -------------------------------------------------------------------------
% Parameter estimation is performed using the custom optimisation wrapper:
%
%   M = aFitDCM(DCM)
%
% followed by:
%
%   M.aloglikVLthermGC
%
% The aFitDCM object builds a reduced parameter representation from the DCM
% priors and calls a Variational Laplace routine in generalised coordinates.
% This routine updates posterior means and covariances by comparing observed
% and predicted spectra under a local Gaussian approximation to the
% posterior.
%
% In practice, the optimisation procedure:
%
%   - proposes parameter updates in reduced space
%   - reconstructs the full parameter structure
%   - evaluates the multi-compartment transfer-function prediction
%   - computes spectral prediction error
%   - updates the posterior mean and covariance using local curvature
%   - iterates until a satisfactory fit is reached or a maximum number of
%     refinement attempts has been made
%
% The main outputs are:
%
%   Ep  - posterior parameter estimates
%   CP  - posterior covariance
%   F   - free energy / model evidence proxy
%
%
% WHY USE THIS MODEL?
% -------------------------------------------------------------------------
% This script is intended for situations where a standard single-compartment
% thalamo-cortical neural mass model may be too coarse. The present
% formulation is useful when one wishes to investigate:
%
%   - dendritic versus somatic contributions to spectral features
%   - compartment-specific routing of excitation and inhibition
%   - the role of presynaptic mechanisms in shaping oscillatory responses
%   - richer synaptic dynamics than static conductance changes alone
%
% In this sense, the model is a more biophysically detailed extension of the
% standard TCM, while still retaining a DCM-style inversion framework.
%
%
% OUTPUTS
% -------------------------------------------------------------------------
% For each dataset, the script saves a fitted DCM structure containing:
%
%   DCM.Ep      posterior parameter estimates
%   DCM.Cp      posterior covariance
%   DCM.F       free energy / model evidence proxy
%   DCM.pred    predicted spectra / cross-spectra
%   DCM.w       frequency vector
%   DCM.G       transfer-function / auxiliary output
%   DCM.series  optional reconstructed series and state diagnostics
%
% These outputs can be used for:
%   - subject-level model inversion
%   - posterior predictive checks
%   - mechanistic interpretation of fitted oscillatory spectra
%   - subsequent group-level analyses
%
%
% DEPENDENCIES
% -------------------------------------------------------------------------
% This script requires:
%
%   atcm   - thalamo-cortical modelling package
%   aoptim - optimisation / inference package
%   SPM    - for core DCM/SPM helper functions
%
% Example dependency:
%   atcm: https://github.com/alexandershaw4/atcm
%
%
% NOTES
% -------------------------------------------------------------------------
% - This example is configured as a one-node model for illustration.
% - The model uses a larger hidden-state space than the standard TCM and may
%   therefore require more careful prior specification and fixed-point
%   stabilisation.
% - The transfer-function formulation makes inversion much faster than full
%   brute-force simulation for spectral fitting.
% - The observation model can be tailored to emphasise different compartmental
%   contributions to the measured signal.
% - Because the model is more flexible, it is especially important to check
%   posterior plausibility, parameter identifiability, and fit quality.
%
%
% INPUT
% -------------------------------------------------------------------------
% i : index (or vector of indices) specifying which dataset(s) to fit from
%     the dataset list.
%
%
% AS2020-2025
% Alex Shaw

% EXAMPLE ONE NODE SETUP:
%==========================================================================

% Data & Design
%--------------------------------------------------------------------------
Data.Datasets     = 'TGB.txt';%'MeanSZDatasets.txt';%'AllSZNoMerge.txt'; % textfile list of LFP SPM datasets (.txt)
Data.Design.X     = [];                % design matrix
Data.Design.name  = {'undefined'};     % condition names
Data.Design.tCode = [1];               % condition codes in SPM
Data.Design.Ic    = [1];               % channel indices
Data.Design.Sname = {'PBVE'};            % channel (node) names
Data.Prefix       = 'VLGC_MCTCM_';      % outputted DCM prefix
Data.Datasets     = atcm.fun.ReadDatasets(Data.Datasets);

% Model space - T = ns x ns, where 1 = Fwd, 2 = Bkw
%--------------------------------------------------------------------------
T = [... % this is a 1-node model; nothing to put here...
    0];
F = (T==1);
B = (T==2);
C = [1]';          % input(s)
L = sparse(1,1);

[p]=fileparts(which('atcm.integrate_1'));p=strrep(p,'+atcm','');addpath(p);


% Set up, over subjects
%--------------------------------------------------------------------------
for i = i;%1:length(Data.Datasets)
    
    % Data Naming & Design Matrix
    %----------------------------------------------------------------------
    DCM          = [];
    [fp fn fe]   = fileparts(Data.Datasets{i});
    DCM.name     = [Data.Prefix fn fe];
    
    DCM.xY.Dfile = Data.Datasets{i};  % original spm datafile
    Ns           = length(F);         % number of regions / modes
    DCM.xU.X     = Data.Design.X;     % design matrix
    DCM.xU.name  = Data.Design.name;  % condition names
    tCode        = Data.Design.tCode; % condition index (in SPM)
    DCM.xY.Ic    = Data.Design.Ic;    % channel indices
    DCM.Sname    = Data.Design.Sname; % channel names
    
    
    if exist(DCM.name);
        fprintf('Skipping model %d/%d - already exists!\n( %s )\n',i,length(Data.Datasets),DCM.name);
        continue;
    end
    
    % Extrinsic Connectivity - Model Space
    %----------------------------------------------------------------------
    DCM.A{1} = F;
    DCM.A{2} = B;
    DCM.A{3} = L;
    DCM.B{1} = DCM.A{1} | DCM.A{2};
    DCM.B(2:length(DCM.xU.X)) = DCM.B;
    DCM.C    = C;
    
    % Function Handles
    %----------------------------------------------------------------------
    %DCM.M.f  = @atcm.tc_hilge2;               % model function handle
    DCM.M.f = @atcm.tc_twocmp_stp;
    DCM.M.IS = @atcm.fun.Alex_LaplaceTFwDNew;            % Alex integrator/transfer function
    %DCM.M.IS = @atcm.fun.Alex_LaplaceTFwDNew;
    DCM.options.SpecFun = @atcm.fun.Afft;    % fft function for IS
    
    % Print Progress
    %----------------------------------------------------------------------
    fprintf('Running Dataset %d / %d\n',i,length(Data.Datasets));
    
    % Frequency range of interest
    fq =  [1 90];
    
    % Prepare Data
    %----------------------------------------------------------------------
    DCM.M.U            = sparse(diag(ones(Ns,1)));  %... ignore [modes]
    DCM.options.trials = tCode;                     %... trial code [GroupDataLocs]
    DCM.options.Tdcm   = [300 1300];                   %... peristimulus time
    DCM.options.Fdcm   = fq;                    %... frequency window
    DCM.options.D      = 1;                         %... downsample
    DCM.options.han    = 1;                         %... apply hanning window
    DCM.options.h      = 4;                         %... number of confounds (DCT)
    DCM.options.DoData = 1;                         %... leave on [custom]
    %DCM.options.baseTdcm   = [-200 0];             %... baseline times [new!]
    DCM.options.Fltdcm = fq;                    %... bp filter [new!]
    DCM.options.UseButterband = fq;

    DCM.options.analysis      = 'CSD';              %... analyse type
    DCM.xY.modality           = 'LFP';              %... ECD or LFP data? [LFP]
    DCM.options.spatial       = 'LFP';              %... spatial model [LFP]
    DCM.options.model         = 'tc6';              %... neural model
    DCM.options.Nmodes        = length(DCM.M.U);    %... number of modes
    
    DCM.options.UseWelch      = 1010;
    DCM.options.FFTSmooth     = 0;
    DCM.options.BeRobust      = 0;
    DCM.options.FrequencyStep = 1/2;
    
    DCM.xY.name = DCM.Sname;
    DCM = atcm.fun.prepcsd(DCM);
    DCM.options.DATA = 1 ;

    DCM.xY.y{:}  = agauss_smooth(abs(DCM.xY.y{:}),.6)';
        
    % Subfunctions and default priors
    %----------------------------------------------------------------------
    DCM = atcm.parameters(DCM,Ns);
            
    % other model options
    %----------------------------------------------------------------------
    DCM.M.solvefixed=0;      % 
    DCM.M.x = zeros(1,8,7);  % init state space: ns x np x nstates
    DCM.M.x(:,:,1)=-70;      % init pop membrane pot [mV]

    % --- Expand state vector: 1×8×7  ->  1×8×10
    % Indices (for clarity)
    iVs = 1; iGE = 2; iGI = 3; iGN = 4; iGB = 5; iGM = 6; iGH = 7;
    iVd = 8; iR  = 9; iU  = 10;

    % If M.x existed with 7 states, expand to 10 (Vd, R, uSTP)
    if size(DCM.M.x,3) < 10
        Xnew = zeros(size(DCM.M.x,1), size(DCM.M.x,2), 10);
        Xnew(:,:,1:size(DCM.M.x,3)) = DCM.M.x;

        % Initialise new states
        % Vd starts near Vs; R~0.9 resources available; uSTP ~ baseline use
        Xnew(:,:,iVd) = Xnew(:,:,iVs);       % dendrite voltage ~ soma
        Xnew(:,:,iR)  = 0.9;                 % resources near 1
        Xnew(:,:,iU)  = 0.2;                 % 'use' (will be overridden by U0 prior below)

        DCM.M.x = Xnew;
    end

        
    load([p '/newpoints3.mat'],'pE','pC')

    pE = spm_unvec(spm_vec(pE)*0,pE);

    pC.ID = pC.ID * 0;
    pC.T  = pC.T *0;
    
    pE.J = DCM.M.x(:)*0-1000;    
    %pE.J(1:8) = log([.6 .8 .4 .6 .4 .6 .4 .4]);
    pE.J(1:8) = log([.2 .99 .1 .8 .1 .2 .05 .1]);
    pC.J = pE.J*0;
    %pC.ID = pC.ID + 1/8;
    pE.L = 0;
    pC.a = pC.a*0;

    pE.Gb = pE.H;
    pC.Gb = [1   0   0   0   0   0   0   0;
             0   1   1   0   0   0   0   0;
             0   0   1   0   0   0   0   0;
             0   0   0   1   1   0   0   0;
             0   0   0   0   1   0   0   0;
             0   0   0   0   1   1   0   0;
             0   0   0   0   0   0   0   0;
             0   0   0   0   0   0   1   0]/64;

    %pC.J(1:8)=1/8;
    pC.d(1) = 1/8;
    pC.d(3) = 1/8;

    logit = @(p) log(p./(1-p));

    % === STP (per presyn population: 8×1 each, used across its outputs)
    if ~isfield(pE,'prel'),  pE.prel  = log(0.6)*ones(8,1); end   % baseline release prob (log-space)
    if ~isfield(pE,'tauR'),  pE.tauR  = 0.6*ones(8,1);      end   % R recovery time (s)
    if ~isfield(pE,'tauU'),  pE.tauU  = 0.2*ones(8,1);      end   % u facilitation time (s)
    if ~isfield(pE,'U0'),    pE.U0    = logit(0.2)*ones(8,1); end % baseline 'use' (stored in logit; code does logistic)

    if ~isfield(pC,'prel'),  pC.prel  = ones(8,1)/8; end   % baseline release prob (log-space)
    if ~isfield(pC,'tauR'),  pC.tauR  = ones(8,1)/8;      end   % R recovery time (s)
    if ~isfield(pC,'tauU'),  pC.tauU  = ones(8,1)/8;      end   % u facilitation time (s)
    if ~isfield(pC,'U0'),    pC.U0    = ones(8,1)/8; end % baseline 'use' (stored in logit; code does logistic)


    % === Two-compartment routing & coupling
    if ~isfield(pE,'gc'),      pE.gc      = log(3);       end     % axial coupling (log-space; ~3 mS default)
    if ~isfield(pE,'w_dend'),  pE.w_dend  = [0.8 0.9];   end     % [AMPA_to_dend, NMDA_to_dend] in [0,1]
    if ~isfield(pE,'w_soma'),  pE.w_soma  = [0.9 0.9];   end     % [GABAa_to_soma, GABAb_to_soma] in [0,1]
    
    if ~isfield(pC,'gc'),      pC.gc      = 1/8;       end     % axial coupling (log-space; ~3 mS default)
    if ~isfield(pC,'w_dend'),  pC.w_dend  = [1 1]/8;   end     % [AMPA_to_dend, NMDA_to_dend] in [0,1]
    if ~isfield(pC,'w_soma'),  pC.w_soma  = [1 1]/2;   end     % [GABAa_to_soma, GABAb_to_soma] in [0,1]

    pE.scale = zeros(4,1);
    pC.scale = ones(4,1)/8;

    % PC0.prel   = (0.3^2) * ones(8,1);   % log-space SD ~0.3
    % PC0.tauR   = (0.2^2) * ones(8,1);
    % PC0.tauU   = (0.1^2) * ones(8,1);
    % PC0.U0     = (0.5^2) * ones(8,1);   % logit-space SD
    % PC0.gc     = (0.3^2);               % log-space SD
    % PC0.w_dend = (0.2^2) * ones(1,2);   % direct [0,1] (clamped in model)
    % PC0.w_soma = (0.2^2) * ones(1,2);

    [J, Wsrc]    = atcm.build_twocmp_observer_J(DCM.M, DCM.M.pE);
    pE.J   = log(J);
    pE.J(isinf(pE.J))=-1000;

    pC.J   = J*0; pC.J(find(J)) = 1/8;


    % Make changes here;
    %-----------------------------------------------------------
   
    DCM.M.pE = pE;
    DCM.M.pC = pC;

    % Optimise using AO.m -- a Newton scheme with add-ons and multiple
    % objective functions built in, including free energy
    %----------------------------------------------------------------------
    w   = DCM.xY.Hz;
    Y   = DCM.xY.y{:};
    DCM.M.y  = DCM.xY.y;
    DCM.M.Hz = DCM.xY.Hz;

    ppE = DCM.M.pE;
    ppC = DCM.M.pC;

    fprintf('--------------- STATE ESTIMATION ---------------\n');
    fprintf('Search for a stable fixed point\n');

    DCM.M.x(:,:,iU) = repmat( 1./(1+exp(-DCM.M.pE.U0(:)))', size(DCM.M.x,1), 1 );  % logistic(U0)

    %xx = load([p '/newx.mat']); 
    %DCM.M.x = spm_unvec(xx.x,DCM.M.x);
    %load('init_14dec','x');
    %DCM.M.x = spm_unvec(x,DCM.M.x);

    %x = atcm.fun.alexfixed(DCM.M.pE,DCM.M,1e-10);
    x = atcm.fun.find_fixed_point_robust(DCM.M.pE,DCM.M,1e-10);
    DCM.M.x = spm_unvec(x,DCM.M.x);

    norm(DCM.M.f(DCM.M.x,0,DCM.M.pE,DCM.M))

    fprintf('Finished...\n');
    
          
    fprintf('--------------- PARAM ESTIMATION ---------------\n');
    %fprintf('iteration %d\n',j);

    % Fit with DCM VB routine:
    %----------------------------------------------------------------------
    %[Qp,Cp,Eh,F] = spm_nlsi_GN(DCM.M,DCM.xU,DCM.xY);
    

    % Fit with Variational Laplace in generalised coordinates
    %----------------------------------------------------------------------
    M = aFitDCM(DCM)

    M.aloglikVLthermGC;%([],0);
    M.update_parameters(M.Ep)

    [y,w,G,s] = feval(DCM.M.IS,spm_unvec(M.Ep,DCM.M.pE),DCM.M,DCM.xU);

    numit = 0;
    while cdist(DCM.xY.y{:}',y{:}') > (1/2) && numit < 8
        numit = numit + 1;
        %M.aloglik;
        M.aloglikVLthermGC;%([],0);
        %M.aloglikFE;
        M.update_parameters(M.Ep);
        [y,w,G,s] = feval(DCM.M.IS,spm_unvec(M.Ep,DCM.M.pE),DCM.M,DCM.xU);
    end
    
    Qp = spm_unvec(M.Ep,DCM.M.pE);
    Cp = M.CP;
    F  = M.F;


    % Fit with LM (Free energy estimation):
    %----------------------------------------------------------------------
    % M = aFitDCM(DCM)
    % 
    % M.aloglikFE
    % M.update_parameters(M.Ep)
    % M.aloglikFE
    % M.update_parameters(M.Ep)
    % M.aloglikFE
    %Qp = spm_unvec(M.Ep,DCM.M.pE);
    %Cp = M.CP;


    % save in DCM structures after optim 
    %----------------------------------------------------------------------
    DCM.M.pE = ppE;
    DCM.Ep = Qp;%spm_unvec(M.Ep,DCM.M.pE);
    DCM.Cp = Cp;

    DCM.M.sim.dt  = 1./600;
    DCM.M.sim.pst = 1000*((0:DCM.M.sim.dt:(2)-DCM.M.sim.dt)');

    [y,w,G,s] = feval(DCM.M.IS,DCM.Ep,DCM.M,DCM.xU);

    DCM.pred = y;
    DCM.w = w;
    DCM.G = G;
    DCM.series = s;
    
    %DCM.Cp = atcm.fun.reembedreducedcovariancematrix(DCM,M.CP);
    %DCM.Cp = makeposdef(DCM.Cp);
    DCM.F  = F;%M.FreeEnergyF;
    %DCM.Cp = M.CP;
    save(DCM.name); close all; clear global;
    
end

end
