function RunTCM_Script(i)
% Top level script showing how to apply the thalamo-cortical neural mass
% model decribed in Shaw et al 2020 NeuroImage, to M/EEG data.
%
% Requires atcm (thalamo cortical modelling package) and aoptim
% (optimisation package)
%
% atcm: https://github.com/alexandershaw4/atcm
% aoptim: https://github.com/alexandershaw4/aoptim
%
% Overview of contents
%--------------------------------------------------
% - atcm. contains:
%         - the equations of motion for an 8 pop thalamo-cortical
%           model described by parameterised Morris-Lecar/Hodgkin-Hux
%           conductance ODEs
%         - numerical integration (Euler, RK, Newton-Cotes +) and spectral
%           response functions
%         - lots of helper functions (integration, differentiation,
%           continuation, decomposition methods etc)
% - aoptim contains:
%          - a second order gradient descent optimisation routine that
%          includes both free energy and orther objective functions
%          - n-th order numerical differentiation functions (parallelised)
%
%
% AS2020/21/22 {alexandershaw4[@]gmail.com}

% EXAMPLE ONE NODE SETUP:
%==========================================================================

% Data & Design
%--------------------------------------------------------------------------
Data.Datasets     = 'MeanSZDatasets.txt'; % textfile list of LFP SPM datasets (.txt)
%Data.Datasets = 'GroupMeans.txt';
Data.Design.X     = [];                % design matrix
Data.Design.name  = {'undefined'};     % condition names
Data.Design.tCode = [1];               % condition codes in SPM
Data.Design.Ic    = [1];               % channel indices
Data.Design.Sname = {'V1'};            % channel (node) names
Data.Prefix       = 'TCM_';      % outputted DCM prefix
Data.Datasets     = atcm.fun.ReadDatasets(Data.Datasets);

% For controls only:
%Data.Datasets=Data.Datasets(contains(Data.Datasets,'control'));

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
    DCM.M.f  = @atcm.tc_hilge;               % model function handle
    DCM.M.IS = @atcm.integrate_1;            % Alex integrator/transfer function
    DCM.options.SpecFun = @atcm.fun.Afft;    % fft function for IS
    
    % Print Progress
    %----------------------------------------------------------------------
    fprintf('Running Dataset %d / %d\n',i,length(Data.Datasets));
    
    % Frequency range of interest
    fq = [3 100];
    
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
    
    DCM.options.UseWelch      = 1010;
    DCM.options.FFTSmooth     = 2;
    DCM.options.BeRobust      = 0;
    DCM.options.FrequencyStep = 1;
    
    DCM.xY.name = DCM.Sname;
    DCM = atcm.fun.prepcsd(DCM);
    DCM.options.DATA = 1 ;
    
    % Do a FOOOF and remove power law
    %fprintf('Removing power law from data spectrum\n');
    %m = atcm.fun.c_oof(DCM.xY.Hz,DCM.xY.y{:});
    %DCM.xY.y{:} =  DCM.xY.y{:} - m;
        
    % Subfunctions and default priors
    %----------------------------------------------------------------------
    DCM = atcm.parameters(DCM,Ns);
    
    % if using AOPTIM for inversion, invoke the linear model g(x) output by
    % placing data (DCM.xY.y) in model struct - DCM.M.y
    DCM.M.y  = DCM.xY.y;
    DCM.M.Hz = DCM.xY.Hz;
    
    % If using DCM inversion, select whether to block graph or not
    DCM.M.nograph = 0;
    
    % Feature function for the integrator [NOT USED]
    %----------------------------------------------------------------------
    DCM = atcm.complete(DCM);
    DCM.M.FS = @(x) x(:).^2.*(1:length(x))'.^2;
    imscale = sum(spm_vec(abs(real(DCM.xY.y{:})))) ./ sum(spm_vec(abs(imag(DCM.xY.y{:}))));
    DCM.M.FS = @(x) [real(x) ; imscale*imag(x) ];
    
    % other model options
    %----------------------------------------------------------------------
    DCM.M.solvefixed=0;      % oscillations == no fixed point search
    DCM.M.x = zeros(1,8,7);  % init state space: ns x np x nstates
    DCM.M.x(:,:,1)=-70;      % init pop membrane pot [mV]
    
    % Set Q - a precision operator, increasing with frequency
    %----------------------------------------------------------------------
    y  = spm_vec(DCM.xY.y{1});
    w  = spm_vec(DCM.xY.Hz);
    Qw = diag(DCM.xY.y{:}./max(DCM.xY.y{:}));
    Qa = Qw;
    Nf = length(w);
    Q  = {spm_Q(1/2,Nf,1)*diag(DCM.M.Hz)*spm_Q(1/2,Nf,1)};
    Qw = Qw * Q{:};
    
    % simulation / integration parameters
    %----------------------------------------------------------------------
    DCM.M.sim.dt  = 1./300;
    DCM.M.sim.pst = 1000*((0:DCM.M.sim.dt:(1)-DCM.M.sim.dt)');
    DCM.M.burnin  = 0;
    DCM.M.intmethod = 0;
    
    % Input is an ERP
    DCM.M.InputType = 2;
    DCM.M.pE.C = log(.01);
    
    % only interested in real psd rn
    %----------------------------------------------------------------------
    DCM.xY.y{1} = real(DCM.xY.y{1});
    DCM.M.y     = DCM.xY.y;
        
    % Parameters
    %if contains(DCM.name,'MeanSch')
        %fprintf('Patient dataset...\n');
        %P = load('aTCM_MeanSch_MeanDataset');
        load('+atcm/PriorSettings','c');
        
        %MP = spm_unvec( (spm_vec(P.M.Ep)+spm_vec(C.M.Ep))/2, P.M.Ep);
        
        %DCM.M.pE = spm_unvec(MP,P.DCM.M.pE);
        DCM.M.pE = c.pE;
        DCM.M.pC = c.pC;
        DCM.M.x  = c.x;
        
        DCM.M.pE.L = 0;
        
        DCM.M.DoHamming = 1;
        DCM.M.pC.C=1/8;
        
        %DCM.M.pC.Gsc([1 2 4 8])=1/8;
        %DCM.M.pC.H(8,8)=1/8;
        %DCM.M.pC.H(6,6)=1/8;
        %DCM.M.pC.H(5,5)=1/8;
        %DCM.M.pC.H(1,1)=1/8;
        
        %DCM.M.pC.S = ones(1,8)/8;
        %DCM.M.pC.Gsc = [1 1 0 1 0 0 0 1]/8;
        
    %else 
    %    fprintf('Control dataset...\n');
    %    P = load('aTCM_MeanCon_MeanDataset');
    %    DCM.M.pE = spm_unvec(P.M.Ep,P.DCM.M.pE);
    %    DCM.M.pC = P.DCM.M.pC;
    %    DCM.M.x  = P.DCM.M.x;
    %end
    
    %EX = load('NEWP','Ep')
    %DCM.M.pE = EX.Ep;
    
    x=load('+atcm/Nov31')
    DCM.M.pE = x.Ep;
    DCM.M.pE.L=-2;
    
    % Optimise using AO.m
    %----------------------------------------------------------------------

    ppE = DCM.M.pE;
    
    M = AODCM(DCM);
    
    % Bias and feature selection - ensuring FS(y) remains smooth
    M.opts.Q  = spm_Q(1/2,Nf,1)*diag(DCM.M.Hz)*spm_Q(1/2,Nf,1);
    Q = atcm.fun.QtoGauss(DCM.xY.y{:},2);
    Q = 1+rescale(Q.*DCM.xY.Hz);
    M.opts.Q = Q;
    
    %M.opts.FS = @(x) [(1+atcm.fun.makef(w(:),50,2,8)).*real(sqrt(denan(x))); denan(std(x)./mean(x)) ];
    
    M.opts.FS = @(x) [real(sqrt(denan(x))); denan(std(x)./mean(x)) ];
    
    % opt set 1.
    M.opts.EnforcePriorProb=0; % forcibly constrain parameters to within prior dist
    M.opts.ismimo      = 0;        % compute dfdp elementwise on vector-output function
    M.opts.doparallel  = 1;    % use parfor loops when poss, incl for df/dx
    M.opts.hyperparams = 0;   % hyperparameter tuning
    M.opts.fsd         = 0;         % fixed-step for derivatives
    M.opts.corrweight  = 0;  % weight log evidence by correlation
    M.opts.inner_loop  = 2;
    
    M.opts.objective           = 'gauss';%'rmse';%'rmse';%'fe';%'mvgkl';%'jsd';%'mvgkl';%'qrmse'; % objective (error) function
    M.opts.criterion           = -inf;%-1000;%1e-3;
    M.opts.isGaussNewton       = 0;
    M.opts.factorise_gradients = 1;
    M.opts.normalise_gradients = 0;
    
    M.opts.hypertune       = 1;
    M.opts.memory_optimise = 1;
    M.opts.rungekutta      = 6;
    M.opts.updateQ         = 1; % do a grd ascent on Q but also weight by residual
    M.opts.crit            = [0 0 0 0];
    M.opts.do_gpr          = 0;
    
    %M.opts.variance_estimation = 0;
    M.opts.ismimo = 1;
    %M.opts.isGaussNewton = 1;
    %M.opts.DoMAP_Bayes = 1;
    %M.opts.docompare=1;
    
    M.opts.userplotfun = @aodcmplotfun;
    M.opts.isGaussNewtonReg=1;
    %M.opts.WeightByProbability=1;
    
    M.default_optimise([7],[44])
                 
    % reinstate the actual priors before saving
    DCM.M.pE = ppE;
    
    DCM.Ep = spm_unvec(M.Ep,DCM.M.pE);
    DCM.Cp = atcm.fun.reembedreducedcovariancematrix(DCM,M.CP);
    DCM.Cp = makeposdef(DCM.Cp);
    DCM.F = M.F;
    
    save(DCM.name); close; clear global;
    
    close all; drawnow;
end