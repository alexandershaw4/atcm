
% EXAMPLE TWO NODE SETUP:
%==========================================================================

% Data & Design
%--------------------------------------------------------------------------
Data.Datasets     = 'KETDEP.txt';  % textfile list of LFP SPM datasets (.txt)
Data.Datasets     = atcm.fun.ReadDatasets(Data.Datasets);
Data.Design.X     = [];              % std/dev
Data.Design.name  = {'Undefined'};         % condition names
Data.Design.tCode = [1];             % condition codes in SPM
Data.Design.Ic    = [1 2];          % channel indices
Data.Design.Sname = {'Frontal' 'Parietal'};  % channel (node) names
Data.Prefix       = 'mTCM_';      % outputted DCM prefix

% Model space - T = ns x ns, where 1 = Fwd, 2 = Bkw
%--------------------------------------------------------------------------
T = [...
    0 1 ;
    2 0 ];

F = (T==1);
B = (T==2);
C = [1 1]';      % inputs
L = sparse(2,2); 


% THE REST IS COMMON ACROSS SETUPS: DONT EDIT, APART FROM 'PREPARE DATA'
%--------------------------------------------------------------------------
modelstore = fileparts(which('atcm.integrate3'));

% Set up, over subjects
%--------------------------------------------------------------------------
for s = 1:length(Data.Datasets)
    
    % Data Naming & Design Matrix
    %----------------------------------------------------------------------
    DCM          = [];
    [fp fn fe]   = fileparts(Data.Datasets{s});
    DCM.name     = [Data.Prefix fn fe];
    
    DCM.xY.Dfile = Data.Datasets{s};  % original spm datafile
    Ns           = length(F);         % number of regions / modes
    DCM.xU.X     = Data.Design.X;     % design matrix
    DCM.xU.name  = Data.Design.name;  % condition names
    tCode        = Data.Design.tCode; % condition index (in SPM)
    DCM.xY.Ic    = Data.Design.Ic;    % channel indices
    DCM.Sname    = Data.Design.Sname; % channel names
        
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
    DCM.M.f  = @atcm.tc_hilge;               % model function handle
    DCM.M.IS = @atcm.integrate3;              % Alex integrator/transfer function
    DCM.options.SpecFun = @atcm.fun.Afft;  % fft function for IS
    
    % Print Progress
    %----------------------------------------------------------------------
    fprintf('Running Dataset %d / %d\n',s,length(Data.Datasets));
    
    fq = [4 90];
    
    % Prepare Data
    %----------------------------------------------------------------------
    DCM.M.U            = sparse(diag(ones(Ns,1)));  %... ignore [modes]
    DCM.options.trials = tCode;                     %... trial code [GroupDataLocs]
    DCM.options.Tdcm   = [300 1300];                   %... peristimulus time
    DCM.options.Fdcm   = fq;                    %... frequency window
    DCM.options.D      = 1;                         %... downsample
    DCM.options.han    = 0;                         %... apply hanning window
    DCM.options.h      = 4;                         %... number of confounds (DCT)
    DCM.options.DoData = 1;                         %... leave on [custom]
    %DCM.options.baseTdcm   = [-200 0];             %... baseline times [new!]
    DCM.options.Fltdcm = fq;                    %... bp filter [new!]

    DCM.options.analysis      = 'CSD';              %... analyse type
    DCM.xY.modality           = 'LFP';              %... ECD or LFP data? [LFP]
    DCM.options.spatial       = 'LFP';              %... spatial model [LFP]
    DCM.options.model         = 'tc6';              %... neural model
    DCM.options.Nmodes        = length(DCM.M.U);    %... number of modes

    % Alex additions - 1010 = use atcm.fun.AFFT.m
    DCM.options.UseWelch      = 1010;
    DCM.options.FFTSmooth     = 4;
    DCM.options.UseButterband = fq;
    DCM.options.BeRobust      = 0;
    DCM.options.FrequencyStep = 1;        % use .5 Hz steps
    
    DCM.xY.name = DCM.Sname;
    DCM = atcm.fun.prepcsd(DCM);
    DCM.options.DATA = 1 ;      
    
    % Subfunctions
    %----------------------------------------------------------------------
    %DCM = atcm.parameters(DCM,Ns,'Priors2021c');       % gets latet priors for tc nmm     
    DCM = atcm.parameters(DCM,Ns,[modelstore '/+fun/Priors2021b']);
    
    DCM.M.pE.gaba = zeros(1,8);
    DCM.M.pC.gaba = zeros(1,8)+1/16;
    
    % if using AOPTIM for inversion, invoke the linear model g(x) output by
    % placing data (DCM.xY.y) in model struct - DCM.M.y
    DCM.M.y  = DCM.xY.y;
    DCM.M.Hz = DCM.xY.Hz;
    
    % If using DCM inversion, select whether to block graph or not
    DCM.M.nograph = 0;
    
    % Final options for integrator
    DCM.M.fmethod = 'none';
    DCM.M.DoEnv   = 0;
    
    DCM.M.ncompe=0;
    DCM.M.envonly=1;
    DCM.M.EnvLFP=1;
    DCM.M.burnin = 300;
    DCM.M.solvefixed=1;
    DCM.M.DoHamming=0;
    DCM.M.LFPsmooth=12;
    DCM.M.usesmoothkernels=0;
    DCM.M.intmethod=2;
    DCM.M.IncDCS=0;
    
    DCM.M.pC.d = zeros(8,1);
    DCM.M.pE.L = repmat(-1.75,[Ns 1]);
    DCM.M.pC.CV = zeros(1,8);
    DCM.M.pC.T = repmat([1 1 1 1]/16,[Ns 1]);
            
    DCM.M.InputType=0; % NOT OSCILLATION
        
    %X = load([modelstore '/+fun/Priors2021a.mat']);
    %DCM.M.pC.H = DCM.M.pC.H + (X.pC.H/2);
    
    % Feature function for the integrator
    %DCM.M.FS = @(x) x(:).^2.*(1:length(x))'.^2;
    DCM = atcm.complete(DCM);
    DCM.M.FS = @(x) x(:).^2.*(1:length(x))'.^2;
    
    % oscillations == no fixed point search
    DCM.M.solvefixed=0;
    DCM.M.x = zeros(Ns,8,7);
    DCM.M.x(:,:,1)=-50;
    DCM.M.ncompe = 47;
    DCM.M.pC.CV = ones(1,8)/8;
    DCM.M.pC.J([2 4])=1/8;
    DCM.M.pC.S = ones(1,8)/16;
        
    %X = load([modelstore '/+fun/NewControlPriors'],'pE');
    %DCM.M.pE = X.pE;
    
    % Optimise BASLEINE                                                  1
    %----------------------------------------------------------------------
    M = AODCM(DCM);
    
    % opt set 1.                   - DONT CHANGE THESE! - 
    M.opts.EnforcePriorProb=0;   % don't force a bayesian approach
    M.opts.ismimo=0;             % don't compute gradients as a MIMO
    M.opts.doparallel=1;         % compute derivatives in parallel (parfor)
    M.opts.hyperparams=1;        % include hyperparameter tuning of precision
    M.opts.fsd=0;                % don't use a fixed step derivative in finite diff computation
    
    %w = DCM.xY.Hz;
    %M.opts.Q=spm_Q(1/2,length(w),1)*diag(w)*spm_Q(1/2,length(w),1);
    
    M.default_optimise([1 3 1],[15 4 4]);
    
    save(DCM.name); close; clear global;    

end