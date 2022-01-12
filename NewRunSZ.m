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
%Data.Datasets     = 'Meansets.txt';
%Data.Datasets     = 'PMP_Datasets.txt';  % textfile list of LFP SPM datasets (.txt)
Data.Design.X     = [];              % std/dev
Data.Design.name  = {'undefined'};         % condition names
Data.Design.tCode = [1];             % condition codes in SPM
Data.Design.Ic    = [1];             % channel indices
Data.Design.Sname = {'V1'};         % channel (node) names
Data.Prefix       = 'nfinalTCM_';      % outputted DCM prefix
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
    DCM.M.f  = @atcm.tc_hilge;               % model function handle
    DCM.M.IS = @atcm.integrate_1;            % Alex integrator/transfer function
    DCM.options.SpecFun = @atcm.fun.Afft;    % fft function for IS
    
    % Print Progress
    %----------------------------------------------------------------------
    fprintf('Running Dataset %d / %d\n',s,length(Data.Datasets));
    
    fq = [3 90];
    
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

    DCM.options.analysis      = 'CSD';              %... analyse type
    DCM.xY.modality           = 'LFP';              %... ECD or LFP data? [LFP]
    DCM.options.spatial       = 'LFP';              %... spatial model [LFP]
    DCM.options.model         = 'tc6';              %... neural model
    DCM.options.Nmodes        = length(DCM.M.U);    %... number of modes

    % Alex additions - 1010 = use atcm.fun.AFFT.m
    DCM.options.UseWelch      = 1010;
    DCM.options.FFTSmooth     = 2;
    %DCM.options.UseButterband = fq;
    DCM.options.BeRobust      = 0;
    DCM.options.FrequencyStep = 1;        % use .5 Hz steps
            
    DCM.xY.name = DCM.Sname;
    DCM = atcm.fun.prepcsd(DCM);
    DCM.options.DATA = 1 ;      
        
    %yx = fit(DCM.xY.Hz.',DCM.xY.y{1},'Gauss4');
    %DCM.xY.y{1} = yx(DCM.xY.Hz);
    
    DCM.xY.yorig = DCM.xY.y;
    
    tmp = real(DCM.xY.y{1});
    w   = DCM.xY.Hz;
    c   = atcm.fun.c_oof(w(:),tmp,'exp2');% constrained exponent
    m   = fit(w(:),real(tmp-c),'Gauss4'); % gmm on residual
    r   = tmp - (m(w) + c);
    tmp = (c/4) + m(w) ;%+ (r/4);
    DCM.xY.y{1} = tmp;
    
    DCM.xY.c = c;
    DCM.xY.r = r;
    DCM.xY.m = m;
    
    
    % Subfunctions and priors
    %----------------------------------------------------------------------
    %DCM = atcm.parameters(DCM,Ns,'Priors2021c');       % gets latet priors for tc nmm     
    DCM = atcm.parameters(DCM,Ns,'~/Dropbox/code/atcm/+atcm/+fun/Priors2021b');
    
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
            
    DCM.M.InputType=2; % NOT OSCILLATION
        
    X = load('~/Dropbox/code/atcm/+atcm/+fun/Priors2021a.mat');
    DCM.M.pC.H = DCM.M.pC.H + (X.pC.H/2);
        
    % Feature function for the integrator
    %----------------------------------------------------------------------
    %DCM.M.FS = @(x) x(:).^2.*(1:length(x))'.^2;
    DCM = atcm.complete(DCM);
    DCM.M.FS = @(x) x(:).^2.*(1:length(x))'.^2;
    
    imscale = sum(spm_vec(abs(real(DCM.xY.y{:})))) ./ sum(spm_vec(abs(imag(DCM.xY.y{:}))));
    DCM.M.FS = @(x) [real(x) ; imscale*imag(x) ];

    % other model options
    %----------------------------------------------------------------------
    DCM.M.solvefixed=0;      % oscillations == no fixed point search
    DCM.M.x = zeros(1,8,7);  % init state space: ns x np x nstates
    DCM.M.x(:,:,1)=-70;      % init pop membrane pot [mV]
    DCM.M.ncompe =0;         % no component analysis
    DCM.M.pC.CV = ones(1,8)/8; 
    DCM.M.pC.J([2 4])=1/8;
    DCM.M.pC.S = ones(1,8)/16;
    
    DCM.M.pE.L = -2.5;
    DCM.M.pC.R = [1 1 1]/8;
    DCM.M.pE.R = [0 0 0];
    
    DCM.M.ppE=DCM.M.pE;
    
    DCM.M.pE.Ly = 2;
    DCM.M.pC.Ly = 1/8;
    
    DCM.M.pE.C = 0;
    DCM.M.pC.C = 1/8;
    
    % Set Q - a precision operator, using central freqs for gauss fit
    %----------------------------------------------------------------------
    y  = spm_vec(DCM.xY.y{1});
    w  = spm_vec(DCM.xY.Hz);
    Qw = diag(DCM.xY.y{:}./max(DCM.xY.y{:}));
    %Qw = 2*exp(-((w-m.b1)/(m.c1/4)).^2) + 2*exp(-((w-m.b2)/(m.c2/4)).^2) + ... 
    %     2*exp(-((w-m.b3)/(m.c3/4)).^2) + 2*exp(-((w-m.b4)/(m.c4/4)).^2);
    %Qw = diag(Qw);
    Nf = length(w);
    Q  = {spm_Q(1/2,Nf,1)*diag(DCM.M.Hz)*spm_Q(1/2,Nf,1)};
    Qw = Qw * Q{:};

    % New 2021: Try fititng with only the synaptioc
    % parameters corresponding to the CMC13 connections!
    %----------------------------------------------------------------------

    pC = spm_unvec(spm_vec(DCM.M.pC)*0,DCM.M.pC);

    pC.H = [1  1  1  0  0  0  0  1;
            1  1  1  0  0  0  0  0;
            1  1  1  0  0  0  0  0;
            1  1  0  1  1  0  0  0;
            0  1  0  1  0  0  0  0;
            0  0  0  1  1  0  0  0;
            0  0  0  0  0  0  0  0;
            0  0  0  0  0  1  1  1]/8;

    pC.Hn= [0  1  0  0  0  0  0  1;
            1  1  0  0  0  0  0  0;
            1  1  0  0  0  0  0  0;
            1  1  0  0  0  0  0  0;
            0  0  0  1  0  0  0  0;
            0  0  0  1  0  0  0  0;
            0  0  0  0  0  0  0  0;
            0  0  0  0  0  1  0  0]/8;

    %DCM.M.pE.T = repmat(DCM.M.pE.T,[8 1]);
    pC.T = [1    0    1    0;
            1    0    1    0;
            0    1    0    1;
            1    0    1    0;
            0    1    0    1;
            1    0    1    0;
            0    1    0    0;
            1    0    1    0]/8;
         
    pC.L = 1/8;
    pC.R = [1 1 1]/8;
    
    DCM.M.pC = pC;
                
    DCM.M.fmethod='none';
    DCM.M.pE.L = -2.5;
        
    DCM.M.InputType=2;
    DCM.M.burnin=300;
        
    % new priors (again...)
    x = load('Priors_22Oct_b','Ep');
    DCM.M.pE = x.Ep;
    DCM.M.pE.L = 0;%-2;
    DCM.M.pE.J([1 4])=-1000;
        
    DCM.M.pE.ID = zeros(1,8);
    DCM.M.pC.ID = ones(1,8)/18;
        
    DCM.M.pC.Ly = 0;
    DCM.M.pE.Ly = -2;
    DCM.M.pE.ID(1)=-0.4;
    
    %DCM.M.pE.f = [0 0];
    %DCM.M.pC.f = [0 0]+1/8;
    
    %DCM.M.pE.b = [0;0];
    %DCM.M.pC.b = [1;1]/8;
        
    DCM.M.y  = DCM.xY.y;

    DCM.M.pE.L   = -.5;
    DCM.M.pC.L   = 4;
    %DCM.M.IncDCS = 1;
    
    %DCM.M.pE.b = [-4;-1];
    %DCM.M.pC.b = [1; 1]/8;
    %DCM.M.pC.b = [0;0];
    %DCM.M.pC.d = ones(8,1)/8;
        
    %DCM.M.pC.E = 1/8;
    %DCM.M.pC.S = ones(1,8)/16;
        
    % keep for saving after
    ppE = DCM.M.pE;
            
    DCM.M.pC   = spm_unvec(spm_vec(DCM.M.pC)./12,DCM.M.pC);
    DCM.M.pC.L = 4;
        
    % Optimise --- 1                                                           1
    %----------------------------------------------------------------------
    M = AODCM(DCM);

    % Bias and feature selection
    M.opts.Q  = real(Qw);  
    M.opts.FS = @(x) [real(sqrt(x))];
    M.opts.FS = @(x) real( spm_vec(atcm.fun.Pf2VMD(x,3)) );
        
    % opt set 1.
    M.opts.EnforcePriorProb=0; % forcibly constrain parameters to within prior dist
    M.opts.ismimo=0;        % compute dfdp elementwise on vector-output function
    M.opts.doparallel=1;    % use parfor loops when poss, incl for df/dx
    M.opts.hyperparams=1;   % hyperparameter tuning
    M.opts.fsd=0;           % fixed-step for derivatives
    M.opts.corrweight = 0;  % weight log evidence by correlation 
        
    M.opts.objective = 'rmse'; % objective (error) function
    M.opts.criterion = 1e-3;
    
    %M.opts.order = 4;
    
    M.default_optimise([7],[28])
    
    % afterward, use AODCM object to loop through the optimisation steps
    % for a visualisation:
    %for i = 1:28; 
    %    dydp(i,:) = spm_vec(M.opts.fun(M.history.p{i}));
    %end
    
    Morig=M;
    
    % Extract fit and run again
    Ep = spm_unvec(M.Ep,DCM.M.pE);
    DCM.M.pE = Ep;

    % Optimise --- 2                                                         1
    %----------------------------------------------------------------------
    M = AODCM(DCM);
    
    % Bias and feature selection
    M.opts.Q  = Qw;  
    M.opts.FS = @(x) real(sqrt(x));
    M.opts.FS = @(x) real( spm_vec(atcm.fun.Pf2VMD(x,3)) );

    % opt set 1.
    M.opts.EnforcePriorProb=0; % forcibly constrain parameters to within prior dist
    M.opts.ismimo=0;        % compute dfdp elementwise on vector-output function
    M.opts.doparallel=1;    % use parfor loops when poss, incl for df/dx
    M.opts.hyperparams=0;   % hyperparameter tuning
    M.opts.fsd=1.2;           % fixed-step for derivatives
    M.opts.corrweight = 0;  % weight log evidence by correlation 
    
    M.opts.objective = 'rmse'; % objective (error) function
    M.opts.criterion = 1e-3;

    M.default_optimise([7],[8])
      
    % reinstate the actual priors before saving
    DCM.M.pE = ppE;
    
    save(DCM.name); close; clear global;    
    
end