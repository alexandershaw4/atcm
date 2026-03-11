function RunTCM_Script_EigenVL_2026(i)
% Top-level pipeline for fitting a conductance-based thalamo-cortical neural
% mass model (TCM) to M/EEG cross-spectral data using an eigendecomposition-
% based transfer-function formulation and Variational Laplace optimisation.
%
% OVERVIEW
% -------------------------------------------------------------------------
% This script implements a Dynamic Causal Modelling (DCM)-style workflow for
% electrophysiological data, specialised to the thalamo-cortical model
% described in Shaw et al. (2020, NeuroImage). In contrast to brute-force
% time-domain integration, this version uses local linearisation of the
% nonlinear neural mass equations around a stable operating point and
% evaluates spectral responses using an eigenmode decomposition of the
% linearised system.
%
% The pipeline is:
%
%   1) Load and prepare M/EEG spectral data
%   2) Specify a thalamo-cortical generative model and connectivity pattern
%   3) Set model priors and initialise the hidden state space
%   4) Search for a stable fixed point of the neural mass equations
%   5) Linearise the system around that operating point
%   6) Eigendecompose the local Jacobian / state transition operator
%   7) Construct a modal transfer function in the frequency domain
%   8) Compare predicted and observed spectra
%   9) Optimise model parameters using Variational Laplace in generalised
%      coordinates
%  10) Save fitted posteriors, predictions, and diagnostics into a DCM struct
%
%
% MODEL
% -------------------------------------------------------------------------
% The underlying generative model is a conductance-based thalamo-cortical
% neural mass model, implemented here with:
%
%   DCM.M.f = @atcm.tc_hilge2
%
% The hidden neuronal states x evolve according to:
%
%   dx/dt = f(x,u,P,M)
%
% where P contains biophysical parameters and u denotes exogenous input. The
% observable signal y is generated from a weighted readout of the hidden
% states after local linearisation around a fixed point. This allows the
% model to retain a mechanistic interpretation while being evaluated
% efficiently in the frequency domain.
%
% In this example the model is a simple one-node thalamo-cortical circuit,
% but the same framework generalises to larger multi-node architectures with
% forward, backward, and lateral connections.
%
%
% DATA PREPARATION
% -------------------------------------------------------------------------
% The script expects a text file containing a list of SPM-format M/EEG
% datasets. These are converted into DCM-compatible cross-spectral form and
% prepared for model inversion.
%
% Preparation includes:
%   - selecting trials / conditions
%   - choosing channels and source labels
%   - restricting the frequency range of interest
%   - estimating and smoothing spectra
%   - assembling the DCM.xY observation structure
%
% Spectral preparation is handled using helper functions such as:
%
%   atcm.fun.prepcsd
%
% so that observed spectra and model predictions share a common
% representation.
%
%
% EIGENMODE TRANSFER-FUNCTION FORMULATION
% -------------------------------------------------------------------------
% The central feature of this version is the use of an eigendecomposition-
% based transfer-function routine:
%
%   DCM.M.IS = @<eigenmode transfer routine>
%
% Rather than evaluating the transfer function solely through direct
% inversion of the resolvent at each frequency, this approach decomposes the
% local linearised system into its constituent dynamical modes.
%
% After fixed-point evaluation, the nonlinear model is locally approximated
% as:
%
%   dx/dt ≈ A x + B u
%
% where A is the Jacobian of the flow with respect to the states and B is
% the input Jacobian. The local dynamics are then decomposed as:
%
%   A = V * Lambda * V^(-1)
%
% where:
%
%   V       - eigenvectors / modal basis
%   Lambda  - eigenvalues / modal poles
%
% This allows the transfer function to be expressed in modal form, so that
% the spectral response becomes a weighted sum of damped dynamical modes.
% Each mode contributes according to its pole location and its projection
% onto the input and observation spaces.
%
% Conceptually, this means that the model prediction can be interpreted not
% only as a spectrum, but as a superposition of underlying oscillatory /
% dynamical motifs, each with its own characteristic frequency and damping.
%
%
% WHY USE EIGENMODE DECOMPOSITION?
% -------------------------------------------------------------------------
% The eigendecomposition transfer formulation is useful because it provides:
%
%   - a mechanistically interpretable decomposition of the spectral response
%   - direct access to modal frequencies and damping rates
%   - an efficient route to spectral prediction after linearisation
%   - insight into which state-space modes dominate specific peaks or bands
%
% In practice, this makes it especially useful for:
%   - interpreting alpha / beta / spindle-like modes
%   - identifying dominant local circuit motifs
%   - linking fitted spectra to eigenvalues of the linearised dynamics
%   - visualising how parameters reshape oscillatory modes
%
% Relative to a plain resolvent evaluation, the eigenmode formulation is
% particularly attractive when the goal is not just fitting spectra but also
% understanding the modal structure of the fitted thalamo-cortical circuit.
%
%
% IMPORTANT NOTE ON DELAYS
% -------------------------------------------------------------------------
% In its simplest form, the eigendecomposition approach diagonalises the
% instantaneous local Jacobian A and therefore provides a clean modal
% description of the undelayed linearised system. If explicit propagation
% delays are included in the transfer formulation, the effective system
% matrix may become frequency-dependent, in which case the modal structure is
% also frequency-dependent.
%
% Accordingly, the precise interpretation of the eigenmodes depends on the
% implementation used:
%
%   - delay-free / instantaneous modes:
%       one global eigendecomposition of A
%
%   - delay-aware modes:
%       frequency-specific eigendecomposition of an effective operator
%
% The present eigendecomposition approach is therefore best understood as a
% modal analysis of the local linearised system, and may serve as either a
% standalone transfer function or an interpretability layer alongside a more
% explicitly delay-aware Laplace-domain transfer function.
%
%
% PRIORS AND PARAMETERISATION
% -------------------------------------------------------------------------
% Priors are defined using:
%
%   DCM = atcm.parameters(DCM,Ns)
%
% and may then be adjusted manually to constrain or emphasise particular
% parameter families. The prior mean structure (pE) and covariance / precision
% structure (pC) govern parameters such as:
%
%   - observation weights
%   - synaptic / conductance gains
%   - intrinsic and extrinsic couplings
%   - spectral smoothing and delay-related terms
%   - channel / source scaling constants
%
% These priors serve both a biological and numerical role: they regularise
% inference and encourage the inversion to remain within plausible dynamical
% regimes.
%
%
% FIXED POINT / OPERATING POINT SEARCH
% -------------------------------------------------------------------------
% Before eigendecomposition, the script estimates a stable operating point
% x* satisfying:
%
%   f(x*,0,P,M) = 0
%
% using a robust fixed-point routine:
%
%   atcm.fun.find_fixed_point_robust
%
% This step is essential because the Jacobian A, and hence the eigenmodes
% themselves, depend on the operating point. In this sense, the modal
% decomposition is local: it characterises the dynamics of the system near
% the inferred baseline state.
%
%
% OPTIMISATION
% -------------------------------------------------------------------------
% Model parameters are optimised using:
%
%   M = aFitDCM(DCM)
%
% followed by:
%
%   M.aloglikVLthermGC
%
% The aFitDCM wrapper converts the structured DCM priors into a reduced
% optimisation vector and calls a Variational Laplace routine in generalised
% coordinates. During inversion, the optimiser:
%
%   - proposes a parameter update
%   - reconstructs the full DCM parameter structure
%   - evaluates the eigenmode transfer-function prediction
%   - compares prediction and data in spectral space
%   - updates the posterior mean and covariance using local curvature
%
% This yields:
%
%   Ep  - posterior parameter means
%   CP  - posterior covariance
%   F   - variational free energy / objective value
%
% The outer fitting loop may repeat the inversion several times if the model
% prediction remains insufficiently close to the observed spectrum, allowing
% iterative refinement of the posterior solution.
%
%
% INTERPRETATION OF FITTED MODES
% -------------------------------------------------------------------------
% A major advantage of this formulation is that the fitted model can be
% interrogated at the level of eigenmodes. For example, one can inspect:
%
%   - eigenvalues in the complex plane
%   - modal damping (real part)
%   - modal frequency (imaginary part)
%   - state participation / eigenvector loadings
%   - contribution of each mode to the predicted spectrum
%
% This makes the approach particularly useful when trying to understand which
% hidden circuit motifs generate observed oscillatory phenomena, and how
% fitted parameters shift the system toward or away from specific dynamical
% regimes.
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
%   DCM.G       transfer-function / modal output
%   DCM.series  optional reconstructed state / signal diagnostics
%
% Depending on the implementation of the eigenmode transfer routine, the
% auxiliary outputs may also include modal quantities such as eigenvalues,
% eigenvectors, residues, or per-mode spectral contributions.
%
%
% DEPENDENCIES
% -------------------------------------------------------------------------
% This script requires:
%
%   atcm   - thalamo-cortical modelling package
%   aoptim - optimisation / inference package
%   SPM    - for DCM/SPM helper functions
%
% Example dependency:
%   atcm: https://github.com/alexandershaw4/atcm
%
%
% NOTES
% -------------------------------------------------------------------------
% - This example uses a simple one-node model for illustration.
% - The eigenmode transfer formulation is best viewed as a local modal
%   approximation around the inferred operating point.
% - It is often more interpretable than a plain transfer-function inversion,
%   because spectral peaks can be related to specific modes.
% - If delays are important, care is needed in interpreting the eigenmodes,
%   as the effective system may become frequency-dependent.
% - Priors, state initialisation, and parameter restrictions may need tuning
%   for different datasets or model variants.
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
Data.Datasets     = 'NewSZ.txt';%'MeanSZDatasets.txt';%'AllSZNoMerge.txt'; % textfile list of LFP SPM datasets (.txt)
Data.Design.X     = [];                % design matrix
Data.Design.name  = {'undefined'};     % condition names
Data.Design.tCode = [1];               % condition codes in SPM
Data.Design.Ic    = [1];               % channel indices
Data.Design.Sname = {'PBVE'};            % channel (node) names
Data.Prefix       = 'EigVL_TCM_';      % outputted DCM prefix
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
    DCM.M.f  = @atcm.tc_hilge2;               % model function handle
    DCM.M.IS = @atcm.fun.Alex_EigenmodeTF   ;            % Alex integrator/transfer function
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
    DCM.options.Tdcm   = [1 2000];                   %... peristimulus time
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
    DCM.options.FrequencyStep = 1;
    
    DCM.xY.name = DCM.Sname;
    DCM = atcm.fun.prepcsd(DCM);
    DCM.options.DATA = 1 ;

    DCM.xY.y{:}  = atcm.fun.agauss_smooth(abs(DCM.xY.y{:}),1)';
        
    % Subfunctions and default priors
    %----------------------------------------------------------------------
    DCM = atcm.parameters(DCM,Ns);
            
    % other model options
    %----------------------------------------------------------------------
    DCM.M.solvefixed=0;      % 
    DCM.M.x = zeros(1,8,7);  % init state space: ns x np x nstates
    DCM.M.x(:,:,1)=-70;      % init pop membrane pot [mV]
        
    load([p '/newpoints3.mat'],'pE','pC')

    pE = spm_unvec(spm_vec(pE)*0,pE);

    pC.ID = pC.ID * 0;
    pC.T  = pC.T *0;
    
    pE.J = pE.J-1000;    
    %pE.J(1:8) = log([.6 .8 .4 .6 .4 .6 .4 .4]);
    pE.J(1:8) = log([.2 .99 .1 .8 .1 .2 .05 .1]);
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

    xx = load([p '/newx.mat']); DCM.M.x = spm_unvec(xx.x,DCM.M.x);
    load('init_14dec','x');
    DCM.M.x = spm_unvec(x,DCM.M.x);

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

    M.aloglikVLthermGC;
    M.update_parameters(M.Ep)

    [y,w,G,s] = feval(DCM.M.IS,spm_unvec(M.Ep,DCM.M.pE),DCM.M,DCM.xU);

    numit = 0;
    while cdist(DCM.xY.y{:}',y{:}') > (1/2) && numit < 8
        numit = numit + 1;
        %M.aloglik;
        M.aloglikVLthermGC;
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
