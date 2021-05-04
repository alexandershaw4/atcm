function NewRunSZ(i)
% Top level script showing how to apply the thalamo-cortical neural mass
% model decribed in Shaw et al 2020 NeuroImage, to M/EEG data. 
%
% Besides the SPM built-ins, all the dependent functions are contained in
% the +atcm/ package.
%
% I include here my script for applying to a dataset from age-matched
% schizophrenia patients and controls.
%
% As of late 2020, I am having more success using an alternative
% optimisation routine 'aoptim', available at https://github.com/alexandershaw4/aoptim
% This routine is similar to the DCM one (spm_nlsi_GN.m) - it still optimises 
% free energy using a gradient descent, but has some extra options such as 
% momentum parameters and line search. 
%
%
% AS2020

% EXAMPLE ONE NODE SETUP:
%==========================================================================
clear global;

% addpath(genpath('~/spm12'));
% addpath(genpath('/home/sapas10/code/atcm/'));
% addpath(genpath('/home/sapas10/code/aoptim/'));

%cd /cubric/scratch/sapas10/tcm/LauSZ/

% Data & Design
%--------------------------------------------------------------------------
Data.Datasets     = 'SZ_Datasets.txt';
%Data.Datasets     = 'PMP_Datasets.txt';  % textfile list of LFP SPM datasets (.txt)
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
    DCM = atcm.parameters(DCM,Ns,'~/code/atcm/+atcm/+fun/Priors2021b');
    
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
    DCM.M.pE.L = -1.75;
    DCM.M.pC.CV = zeros(1,8);
    DCM.M.pC.T = [1 1 1 1]/16;
            
    DCM.M.InputType=1; % OSCILLATION
    
    % load best so far
    %X = load('NewStartPoint.mat','Ep');
    %DCM.M.pE = X.Ep;
    
    X = load('/Users/Alex/code/atcm/+atcm/+fun/Priors2021a.mat');
    DCM.M.pC.H = DCM.M.pC.H + (X.pC.H/2);
    
    % Feature function for the integrator
    %DCM.M.FS = @(x) x(:).^2.*(1:length(x))'.^2;
    DCM = atcm.complete(DCM);
    DCM.M.FS = @(x) x(:).^2.*(1:length(x))'.^2;
    
    % Optimise BASLEINE                                                  1
    %----------------------------------------------------------------------
    M = AODCM(DCM);
    
    % opt set 1.
    M.opts.EnforcePriorProb=1;
    M.opts.ismimo=0;
    M.opts.doparallel=1;
    M.opts.hyperparams=0;
    M.opts.fsd=0;
    
    w = DCM.xY.Hz;
    M.opts.Q=spm_Q(1/2,length(w),1)*diag(w)*spm_Q(1/2,length(w),1);
    
    M.default_optimise([1],[15]);
    
    save(DCM.name); close; clear global;    
    
end