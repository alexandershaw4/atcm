% Top level script showing how to apply the thalamo-cortical neural mass
% model decribed in Shaw et al 2020 NeuroImage, to M/EEG data. 
%
% Besides the SPM built-ins, all the dependent functions are contained in
% the +atcm/ package.
%
% I include here examples of how to set up both single-region and
% multi-region models.
%
% As of late 2020, I am having more success using an alternative
% optimisation routine 'aoptim', available at https://github.com/alexandershaw4/aoptim
% This routine is similar to the DCM one (spm_nlsi)GN.m) - it still optimises 
% free energy using a gradient descent, but has some extra options such as 
% momentum parameters and line search. 
%
%
% AS2020

% EXAMPLE ONE NODE SETUP:
%==========================================================================

% Data & Design
%--------------------------------------------------------------------------
Data.Datasets     = 'SZ_Datasets.txt';  % textfile list of LFP SPM datasets (.txt)
Data.Design.X     = [];              % std/dev
Data.Design.name  = {'undefined'};         % condition names
Data.Design.tCode = [1];             % condition codes in SPM
Data.Design.Ic    = [1];             % channel indices
Data.Design.Sname = {'V1'};         % channel (node) names
Data.Prefix       = 'TCM_';      % outputted DCM prefix
Data.Datasets     = atcm.fun.ReadDatasets(Data.Datasets);

% Model space - T = ns x ns, where 1 = Fwd, 2 = Bkw
%--------------------------------------------------------------------------
T = [... % this is a 1-node model; nothing to put here...
    0];
F = (T==1);
B = (T==2);
C = [1]';          % input(s)
L = sparse(1,1); 


% EXAMPLE EIGHT NODE SETUP:
%==========================================================================

% % Data & Design
% %--------------------------------------------------------------------------
% Data.Datasets     = {'314_27934_VisAAL'};  % textfile list of LFP SPM datasets (.txt)
% Data.Design.X     = [];              % std/dev
% Data.Design.name  = {'Stim'};         % condition names
% Data.Design.tCode = [1];             % condition codes in SPM
% % l / r thal: 77/78
% % l / r calc: 43/44
% % l / r sup occ: 49/50 (dorsal stream)
% % l / r inf occ: 53/54 (central stream)
% Data.Design.Ic    = [77 78 43 44 49 50 53 54];          % channel indices
% Data.Design.Sname = {'LTh' 'RTh' 'LV1' 'RV1' 'LSupOc' 'RSupOv' 'LInfOc' 'RInfOc'};  % channel (node) names
% Data.Prefix       = 'TCM_';      % outputted DCM prefix
% 
% % Model space - T = ns x ns, where 1 = Fwd, 2 = Bkw
% %--------------------------------------------------------------------------
% T = [...
%     0 0 2 0 0 0 0 0;  % th
%     0 0 0 2 0 0 0 0;  % th
%     1 0 0 0 2 0 2 0;  % v1
%     0 1 0 0 0 2 0 2;  % v1
%     0 0 1 0 0 0 0 0;  % supoc
%     0 0 0 1 0 0 0 0;  % supoc 
%     0 0 1 0 0 0 0 0;  % infoc
%     0 0 0 1 0 0 0 0 ];% infoc
% 
% F = (T==1);
% B = (T==2);
% C = [1 1 0 0 0 0 0 0]';      % inputs
% L = sparse(8,8); 


% THE REST IS COMMON ACROSS SETUPS: 
%                                   DONT EDIT, APART FROM 'PREPARE DATA'
%--------------------------------------------------------------------------





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
    DCM.M.f  = @atcm.tc_dev_dev;               % model function handle
    DCM.M.IS = @atcm.integrate3;              % Alex integrator/transfer function
    %DCM.M.IS = 'spm_csd_mtf';              % default DCM transfer function
    DCM.options.SpecFun = @atcm.fun.Afft;  % fft function for IS
    
    % Print Progress
    %----------------------------------------------------------------------
    fprintf('Running Dataset %d / %d\n',s,length(Data.Datasets));
    
    % Prepare Data
    %----------------------------------------------------------------------
    DCM.M.U            = sparse(diag(ones(Ns,1)));  %... ignore [modes]
    DCM.options.trials = tCode;                     %... trial code [GroupDataLocs]
    DCM.options.Tdcm   = [300 1000];                   %... peristimulus time
    DCM.options.Fdcm   = [4 90];                    %... frequency window
    DCM.options.D      = 1;                         %... downsample
    DCM.options.han    = 0;                         %... apply hanning window
    DCM.options.h      = 4;                         %... number of confounds (DCT)
    DCM.options.DoData = 1;                         %... leave on [custom]
    %DCM.options.baseTdcm   = [-200 0];             %... baseline times [new!]
    DCM.options.Fltdcm = [4 90];                    %... bp filter [new!]

    DCM.options.analysis      = 'CSD';              %... analyse type
    DCM.xY.modality           = 'LFP';              %... ECD or LFP data? [LFP]
    DCM.options.spatial       = 'LFP';              %... spatial model [LFP]
    DCM.options.model         = 'tc6';              %... neural model
    DCM.options.Nmodes        = length(DCM.M.U);    %... number of modes

    % Alex additions - 1010 = use atcm.fun.AFFT.m
    DCM.options.UseWelch      = 1010;
    DCM.options.FFTSmooth     = 4;
    DCM.options.UseButterband = [4 90];
    DCM.options.BeRobust      = 0;
    DCM.options.FrequencyStep = 1;        % use .5 Hz steps
    
    DCM.xY.name = DCM.Sname;
    DCM = atcm.fun.prepcsd(DCM);
    DCM.options.DATA = 1 ;      
    
    % Subfunctions
    %----------------------------------------------------------------------
    DCM = atcm.parameters(DCM,Ns);       % gets latet priors for tc nmm     
    
            
    % if using AOPTIM for inversion, invoke the linear model g(x) output by
    % placing data (DCM.xY.y) in model struct - DCM.M.y
    DCM.M.y  = DCM.xY.y;
    DCM.M.Hz = DCM.xY.Hz;
    
    % If using DCM inversion, select whether to block graph or not
    DCM.M.nograph = 0;
    
    % Final options for integrator
    DCM.M.fmethod = 'none';
    DCM.M.DoEnv   = 1;
    
    DCM.M.ncompe=0;
    DCM.M.envonly=1;
    DCM.M.EnvLFP=1;
    DCM.M.IncDCS=0;
    DCM.M.burnin = 300;
    DCM.M.solvefixed=1;
    DCM.M.DoHamming=0;
    DCM.M.LFPsmooth=0;
    DCM.M.usesmoothkernels=0;
    
    % new parameters
    pC = DCM.M.pC;
    V  = spm_unvec(spm_vec(pC)*0,pC);
    
    V.L  = 1/8;
    V.ID = ones(1,8)*0.0156;
    
    V.H([2 3],3)=1/8;
    V.H([1 2],8)=1/8;
    V.H(2,2)=1/8;
    V.H(3,2)=1/8;
    V.H(2,1)=1/8;
    V.H([4 5],5)=1/8;
    V.H(8,6)=1/8;
    V.H(1,1)=1/8;
    V.H(8,8)=1/8;
    
    V.Hn([2 3],2)=1/8;
    V.Hn(3,2)=1/8;
    V.Hn(2,1)=1/8;
    V.Hn(8,6)=1/8;
    V.Hn(1,8)=1/8;
    
    V.H(4,2)=1/8;
    V.H([4 5 6],4)=1/8;
    V.H(4,5)=1/8;    
    V.H([3 5],[5 2])=1/16;
    
    V.H = pC.H;
    V.Hn = pC.Hn;
    
    V.CV = ones(1,8)/8;
    V.S = ones(1,8)/8;
    
    DCM.M.pC=V;
    DCM.M.pE.L=0;
        
    DCM.M.DoHamming=0;
    DCM.M.DoEnv=1;
    DCM.M.LFPsmooth=12;
    DCM.M.usesmoothkernels=0;
    DCM.M.ncompe=0;
    DCM.M.burnin=300;
    
    DCM.M.pE.J([1 4 6 8])=-1000;
    DCM.M.pE.J(4)=log(.6);
    DCM.M.pE.J(1)=log(.4);
    
    DCM.M.InputType=1;
    DCM.M.pE.R(2)=0;
    DCM.M.pC.R(2)=1/16;
    
    %DCM.M.pE.C = [0 0 0 0 0];
    %DCM.M.pC.C = [1 1 1 1 1]/8;
    
    DCM.M.pE.L = -0.25;
    %DCM.M.pE.L = -1;
    
    
    DCM.M.pE.Gsc = zeros(1,8);
    DCM.M.pC.Gsc = ones(1,8)/16;
    
    DCM.M.pC.b = [1;1]/8;
    
    DCM = atcm.complete(DCM);            % complete network specification
    
    % Optimise BASLEINE                                                  1
    %----------------------------------------------------------------------
    M = AODCM(DCM);

    % opt set 1.
    M.opts.EnforcePriorProb=1;
    M.opts.ismimo=0;
    M.opts.doparallel=1;
    M.opts.hyperparams=1;
    M.opts.fsd=0;
    %M.opts.Q=spm_Q(1/2,length(w),1)*diag(w)*spm_Q(1/2,length(w),1);
    %M.opts.FS = @(x) x(:).^2.*(1:length(x))'.^2;  
    
    M.default_optimise([1],[12]);
    
    % set 2
    M.update_parameters(M.Ep);
    M.opts.Q = [];%QQ;
    M.opts.hyperparams=0;
    M.opts.fsd=0;
    w = DCM.xY.Hz;
    % bias both precision (Q) and feature selection (FS) toward gamma
    M.opts.FS =@(x) x(:).^2.*(1:length(x))'.^2;     
    M.opts.Q = spm_Q(1/2,length(w),1)*diag(w)*spm_Q(1/2,length(w),1);
    M.default_optimise([1],[4]);
    
    % set 3
    M.update_parameters(M.Ep);
    M.opts.ismimo=1;
    M.opts.hyperparams=0;
    M.opts.fsd=0;
    M.opts.FS = [];%@(x) x(:).^2.*(1:length(x))'.^2;
    M.default_optimise([1],[4]);
    
    % set 4
    M.update_parameters(M.Ep);
    M.opts.ismimo=0;
    M.opts.hyperparams=1;
    M.opts.fsd=1;
    M.opts.FS = [];
    M.default_optimise([1],[2]);
    
    Ep = spm_unvec(M.Ep,DCM.M.pE);
    save(DCM.name); close; clear global;

    
    
    
end