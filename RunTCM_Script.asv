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
Data.Design.X     = [];                % design matrix
Data.Design.name  = {'undefined'};     % condition names
Data.Design.tCode = [1];               % condition codes in SPM
Data.Design.Ic    = [1];               % channel indices
Data.Design.Sname = {'V1'};            % channel (node) names
Data.Prefix       = 'aTCM_';      % outputted DCM prefix
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
    
    %m = atcm.fun.c_oof(DCM.xY.Hz,DCM.xY.y{:});
    %Y = DCM.xY.y{:} - m;
    %M = fit(DCM.xY.Hz.',Y,'Gauss5');
    %DCM.xY.y{:} = M(DCM.xY.Hz);
        
    % Subfunctions and default priors
    %----------------------------------------------------------------------
    DCM = atcm.parameters(DCM,Ns);
    
    DCM.xY.y{:} = abs(DCM.xY.y{:});
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
    
    % simulation / integration parameters
    %----------------------------------------------------------------------
    DCM.M.sim.dt  = 1./600;
    DCM.M.sim.pst = 1000*((0:DCM.M.sim.dt:(2)-DCM.M.sim.dt)');
    DCM.M.burnin  = 640;
    
    % Input is an ERP
    DCM.M.InputType = 1;

    DCM.M.UseSmooth=1;
    
    % Use a 2-point RK method for integration
    DCM.M.intmethod = 45;

    % No hamming on spectrum
    DCM.M.DoHamming = 0;

    % USE PROVIDED PRIORS!
    x = load('+atcm/TCM_Priors_Latest.mat','pE','pC');
    DCM.M.pE = x.pE;
    DCM.M.pC = x.pC;

    DCM.M.pE.J([1 2 3 4 5 6 7 8]) = log([.6 .8 .4 .6 .4 .6 .2 .2]);
    DCM.M.pE.ID = zeros(1,8);
    DCM.M.pC.ID = ones(1,8)/8;
    DCM.M.pC.Gsc = ones(1,8)/8;
    DCM.M.pC.R = [1 1]/8;
    DCM.M.pE.L=-2.5;

    %DCM.M.pC.S = ones(1,8)/8;
    %DCM.M.pC.J(1:8)=1/8;
    DCM.M.pC.d = DCM.M.pC.d*0;

    % flat priors
    DCM.M.pE = spm_unvec( real(spm_vec(DCM.M.pE)*0), DCM.M.pE);
    DCM.M.pE.J = DCM.M.pE.J-1000;
    DCM.M.pE.J(2)=log(1.1);
    DCM.M.pE.L=-3;

    DCM.M.pE.dd = ones(8,1)*0;
    DCM.M.pC.dd = ones(8,1)/8;
    DCM.M.pE.L = -5;

    % Optimise using AO.m -- a Newton scheme with add-ons and multiple
    % objective functions built in, including free energy
    %----------------------------------------------------------------------
    ppE = DCM.M.pE;
    
    % Construct an AO optimisation object
    M = AODCM(DCM);
    
    % Bias and feature selection - ensuring FS(y) remains smooth
    %Q = (AGenQn(rescale(atcm.fun.makef(w,median(w),2,32),.5,1),16))';
    %Q = AGenQn(hY,8);
    M.opts.Q = full(DCM.xY.Q);
    
    % Feature selection: FS(y)
    %M.opts.RFS = @(x) [real(sqrt(denan(x))); denan(std(x)./mean(x)) ];

    %M.opts.FS = @(x) [real(sqrt(denan(x))); spm_vec(atcm.fun.maxpointsinds(x,length(x))./length(x))];
        
    % Optimisation option set 1.
    M.opts.EnforcePriorProb=0; % forcibly constrain parameters to within prior dist
    M.opts.ismimo      = 1;        % compute dfdp elementwise on vector-output function
    M.opts.doparallel  = 1;    % use parfor loops when poss, incl for df/dx
    M.opts.hyperparams = 1;   % hyperparameter tuning
    M.opts.fsd         = 1;         % fixed-step for derivatives
    M.opts.corrweight  = 0;  % weight log evidence by correlation
    M.opts.inner_loop  = 1;
    
    M.opts.objective           = 'gauss';%'rmse';%'rmse';%'fe';%'mvgkl';%'jsd';%'mvgkl';%'qrmse'; % objective (error) function
    M.opts.criterion           = -inf;%-1000;%1e-3;
    M.opts.isGaussNewton       = 0;
    M.opts.factorise_gradients = 0;
    M.opts.normalise_gradients = 0;
    
    M.opts.hypertune       = 0; % no
    M.opts.memory_optimise = 1;
    M.opts.rungekutta      = 6;
    M.opts.updateQ         = 1; % do a grd ascent on Q but also weight by residual
    M.opts.crit            = [0 0 0 0];
    M.opts.do_gpr          = 0;
    
    M.opts.userplotfun = @aodcmplotfun;
    M.opts.isGaussNewtonReg=0;        % no
    M.opts.order=1;

    %M.opts.orthogradient=1;
    M.opts.orthogradient=1;
        
    M.default_optimise([9],[18])
    
    M.update_parameters(M.Ep);

    M.default_optimise([9],[8])

    % save after first optim loop (because some fail in stage 2)
    %----------------------------------------------------------------------
    DCM.M.pE = ppE;
    DCM.Ep = spm_unvec(M.Ep,DCM.M.pE);
    DCM.Cp = atcm.fun.reembedreducedcovariancematrix(DCM,M.CP);
    DCM.Cp = makeposdef(DCM.Cp);
    DCM.F = M.F;
    save(DCM.name); close; clear global;
    
end
end