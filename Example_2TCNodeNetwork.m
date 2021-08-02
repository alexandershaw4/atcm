function Example_2TCNodeNetwork(i)


% EXAMPLE ONE NODE SETUP:
%==========================================================================
clear global;


% Data & Design
%--------------------------------------------------------------------------
Data.Datasets     = 'nKETDEP.txt';  % textfile list of LFP SPM datasets (.txt)
Data.Datasets     = atcm.fun.ReadDatasets(Data.Datasets);
Data.Design.X     = [0];              % std/dev
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
L = zeros(2); 

modelstore = fileparts(which('atcm.integrate3'));

% Set up, over subjects
%--------------------------------------------------------------------------
for s = i;%1:length(Data.Datasets)
    
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
    DCM.options.Tdcm   = [0 2000];                   %... peristimulus time
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
    DCM.options.model         = 'tc8';              %... neural model
    DCM.options.Nmodes        = length(DCM.M.U);    %... number of modes

    % Alex additions - 1010 = use atcm.fun.AFFT.m
    DCM.options.UseWelch      = 1010;
    DCM.options.FFTSmooth     = 4;
    DCM.options.UseButterband = fq;
    DCM.options.BeRobust      = 0;
    DCM.options.FrequencyStep = 1;        % use .5 Hz steps
    
    % pre = 1:190, bolus = 191:246, infusion = 247:437
    
    % - need to load spm meeg file to get the .info we stored with
    % the original trial definitions - then find the ones that match the
    % window we're interested in:
    
    %[pre]
    D = spm_eeg_load(Data.Datasets{s});
    DCM.options.TrialIndices = D.info(ismember(D.info,1:190));
    
%     % [bolus]:
%     D = spm_eeg_load(Data.Datasets{s});
%     DCM.options.TrialIndices = D.info(ismember(D.info,191:246));
% 
%     % [infusion]:
%     D = spm_eeg_load(Data.Datasets{s});
%     DCM.options.TrialIndices = D.info(ismember(D.info,247:437));

    DCM.xY.name = DCM.Sname;
    DCM = atcm.fun.prepcsd(DCM);
    DCM.options.DATA = 1 ;      
    
    % Get parameters and options
    %----------------------------------------------------------------------
    DCM = atcm.parameters(DCM,Ns,[modelstore '/+fun/Priors2021b']);
        
    % place data (DCM.xY.y) in model struct - DCM.M.y
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
    
    DCM.M.pC.d = zeros(8,1) ;
    DCM.M.pE.L = -1.75;
    DCM.M.pC.CV = zeros(1,8);
    DCM.M.pC.T = repmat([1 1 1 1]/16,[Ns 1]);
            
    DCM.M.pE.gaba = zeros(1,8);
    DCM.M.pC.gaba = ones(1,8)/8;
    
    DCM.M.InputType=2; % NOT OSCILLATION
        
    X = load([modelstore '/+fun/Priors2021a.mat']);
    DCM.M.pC.H = DCM.M.pC.H + repmat(X.pC.H,[1 1 2]);
    
    % Feature function for the integrator
    %----------------------------------------------------------------------
    %DCM.M.FS = @(x) x(:).^2.*(1:length(x))'.^2;
    DCM = atcm.complete(DCM);
    DCM.M.FS = @(x) x(:).^2.*(1:length(x))'.^2;
    
    % append cross correlation function
    tf = @(x) x(:).^2.*(1:length(x))'.^2;;
    %DCM.M.FS = @(x) [(x) ; spm_csd2ccf(tf(x),DCM.xY.Hz,1./32) ];
    imscale = sum(spm_vec(abs(real(DCM.xY.y{:})))) ./ sum(spm_vec(abs(imag(DCM.xY.y{:}))));
    
    % FS = stack real(x) and scaled imag(x)
    DCM.M.FS = @(x) [real(x) ; imscale*imag(x) ];
    
    % Look for a stable point?
    %----------------------------------------------------------------------
    % oscillations == no fixed point search
    DCM.M.solvefixed=0;
    DCM.M.x = zeros(2,8,7);
    DCM.M.x(:,:,1)=-70;
    DCM.M.ncompe =0;
    DCM.M.pC.CV = ones(1,8)/8;
    DCM.M.pC.J([2 4])=1/8;
    DCM.M.pC.S = ones(1,8)/16;
    
    % Final tweaks of parameters and options 
    %----------------------------------------------------------------------
    DCM.M.pE.L = repmat(-2.5,[Ns,1]);
    DCM.M.pC.R = [1 1 1]/8;
    DCM.M.pE.R = [0 0 0];
    DCM.M.pE.L = [3.5 -3.6]';    
    DCM.M.ppE  = DCM.M.pE;
    
    DCM.M.pE.Ly = 2;
    DCM.M.pC.Ly = 1/8;
    
    % Use an exponantial decay curve as noise on the cross-specra only
    % (not autospec)
    DCM.M.pE.c = [1 1; -0.3475 -0.3475];
    DCM.M.pC.c = [1 1; 1 1]/8;
    DCM.M.pE.L(1)=4;
    DCM.M.pE.iL = [0; 0];
    DCM.M.pC.iL = [1; 1];
    DCM.M.DoEnv=1;
        
    % Set Q: precision matrix for the optimiser
    %----------------------------------------------------------------------
    for nn = 1:Ns
        for no = 1:Ns
        
            y  = spm_vec(DCM.xY.y{1}(:,nn,no));
            w  = spm_vec(DCM.xY.Hz);
            [~,LO] = findpeaks(real(smooth(y)),w,'NPeak',4);
            Qw = zeros(size(w))+1;
            for ip = 1:length(LO)
                i0(ip)=atcm.fun.findthenearest(w,LO(ip));
            end
            Qw(i0)=4;
            Qw=diag(Qw);

            QQ{nn,no} = Qw;
        end
    end
    
    Qw = diag([diag(QQ{1,1}); diag(QQ{2,1}); diag(QQ{1,2}); diag(QQ{2,2})]);
    
    % rescale spectra
    DCM.xY.y{1}(:,1,1) = DCM.xY.y{1}(:,1,1)./sum((DCM.xY.y{1}(:,1,1)));
    DCM.xY.y{1}(:,2,2) = DCM.xY.y{1}(:,2,2)./sum((DCM.xY.y{1}(:,2,2)));
    DCM.xY.y{1}(:,1,2) = DCM.xY.y{1}(:,1,1).*conj( DCM.xY.y{1}(:,2,2) );
    DCM.xY.y{1}(:,2,1) = DCM.xY.y{1}(:,2,2).*conj( DCM.xY.y{1}(:,1,1) );
    
    % Optimise [aka invert]                                              1
    %----------------------------------------------------------------------
    M = AODCM(DCM);
    
    M.opts.Q = Qw;
    
    M.opts.FS = DCM.M.FS;
    
    % opt set 1.
    M.opts.EnforcePriorProb=0;
    M.opts.ismimo=0;
    M.opts.doparallel=1;
    M.opts.hyperparams=1;
    M.opts.fsd=0;
    M.opts.corrweight = 1; % weight error by correlation (good for spectra)
    
    % add user-defined plot function
    M.opts.userplotfun = @aodcmplotfun;
    
    %w = DCM.xY.Hz;
    %M.opts.Q=spm_Q(1/2,length(w),1)*diag(w)*spm_Q(1/2,length(w),1);
    
    M.default_optimise([1 3 1],[15 4 4]);
    
    save(DCM.name); close; clear global;    
    
end