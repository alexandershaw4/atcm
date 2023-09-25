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
Data.Datasets     = 'AllSZNoMerge.txt';%'MeanSZDatasets.txt';%'AllSZNoMerge.txt'; % textfile list of LFP SPM datasets (.txt)
Data.Design.X     = [];                % design matrix
Data.Design.name  = {'undefined'};     % condition names
Data.Design.tCode = [1];               % condition codes in SPM
Data.Design.Ic    = [1];               % channel indices
Data.Design.Sname = {'V1'};            % channel (node) names
Data.Prefix       = 'TCMx_';      % outputted DCM prefix
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
    DCM.M.f  = @atcm.tc_hilge2;               % model function handle
    DCM.M.IS = @atcm.integrate_1;            % Alex integrator/transfer function
    DCM.options.SpecFun = @atcm.fun.Afft;    % fft function for IS
    
    % Print Progress
    %----------------------------------------------------------------------
    fprintf('Running Dataset %d / %d\n',i,length(Data.Datasets));
    
    % Frequency range of interest
    fq = [1 90];
    
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
    DCM.options.FFTSmooth     = 1;
    DCM.options.BeRobust      = 0;
    DCM.options.FrequencyStep = 1;
    
    DCM.xY.name = DCM.Sname;
    DCM = atcm.fun.prepcsd(DCM);
    DCM.options.DATA = 1 ;

    % also without the robust fitting to get the residual
    %DCMo = DCM;
    %DCMo.options.BeRobust=0;
    %DCMo = atcm.fun.prepcsd(DCMo);
    %r = DCMo.xY.y{1} - DCM.xY.y{1};
        
    % Subfunctions and default priors
    %----------------------------------------------------------------------
    DCM = atcm.parameters(DCM,Ns);
    
    DCM.xY.y{:} = abs(DCM.xY.y{:});
    
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
        
    % simulation / integration parameters
    %----------------------------------------------------------------------
    DCM.M.sim.dt  = 1./600;
    DCM.M.sim.pst = 1000*((0:DCM.M.sim.dt:(2)-DCM.M.sim.dt)');
    DCM.M.burnin  = 640;
    
    % Input is an oscillation
    DCM.M.InputType = 1;

    % Use a 2-point RK method for integration
    DCM.M.intmethod = 45;

    % No hamming on spectrum
    DCM.M.DoHamming = 0;

    load("newpriors.mat")
    DCM.M.pE = Ep;
    DCM.M.pC = pC;


    DCM.M.pC.S  = ones(1,8)/8;

    DCM.M.pE.T = [0 0 0 0 0 0];
    DCM.M.pC.T = [1 1 1 1 1 1]/8;

    DCM.M.pE.Mh = zeros(1,8);
    DCM.M.pC.Mh = [0 0 0 0 0 1 0 1]*0;%/8;

    DCM.M.pE.Hh = zeros(1,8);
    DCM.M.pC.Hh = [0 0 0 0 0 1 0 1]*0;%/8;

    % remove old params
    rm = {'h' 'm' 'gaba' 'psmooth'};

    for i = 1:length(rm)
        try DCM.M.pE = rmfield(DCM.M.pE,rm{i});end
        try DCM.M.pC = rmfield(DCM.M.pC,rm{i});end
    end
    
    % generate a confounds Q matrix
    Y = DCM.xY.y{:};
    w  = DCM.xY.Hz;

    Q = spm_Q(1/2,length(w),1)*diag(w)*spm_Q(1/2,length(w),1)';    
    Q = Q .* atcm.fun.AGenQn(DCM.xY.y{:},8);
    Q = abs(Q) + AGenQn(diag(Q),8);
    
    % July 2023: try turning everything off except intrinsics and delays:
    pE = DCM.M.pE;
    pC = DCM.M.pC;

    DCM.M.pE   = spm_unvec( real(spm_vec(DCM.M.pE)*0), DCM.M.pE);
    DCM.M.pE.J = pE.J;  
    DCM.M.pE.L = -2;

    DCM.M.pC = spm_unvec(spm_vec(DCM.M.pC)*0,DCM.M.pC);
    DCM.M.pC.H  = pC.H;
    DCM.M.pC.Hn = pC.Hn;
    DCM.M.pC.ID = ones(1,8)*0;%./32;%/16;
    DCM.M.pC.L = 1/8;
    
    DCM.M.pC.H(6,6) = 1/8;
    DCM.M.pC.CT = 1/8;
    DCM.M.pC.TC = 1/8;    

    load('new_priors_6923','pE')
    DCM.M.pE   = pE;
    DCM.M.pE.L = -2;

    DCM.M.pC.S = DCM.M.pC.S + 1/8;
    DCM.M.pC.ID = ones(1,8)/32;
    
    % Optimise using AO.m -- a Newton scheme with add-ons and multiple
    % objective functions built in, including free energy
    %----------------------------------------------------------------------
    w   = DCM.xY.Hz;
    Y   = DCM.xY.y{:};

    DCM.M.y  = DCM.xY.y;
    DCM.M.Hz = DCM.xY.Hz;

    DCM.M.InputType=9;
    DCM.M.pC.d(1:6)=1/8;
    DCM.M.pE.J(1) = log(.4);
    ppE = DCM.M.pE;
    DCM.M.solvefixed=0;
    %DCM.M.x = atcm.fun.solvefixedpoint(DCM.M.pE,DCM.M,[],-70);
            
    % Construct an AO optimisation object
    M = AODCM(DCM);

    M.opts.Q = Q;%diag(w);
    
    % Optimisation option set 1.
    M.opts.EnforcePriorProb    = 0; 
    M.opts.WeightByProbability = 0;

    M.opts.ismimo      = 1;     
    M.opts.doparallel  = 1;    
    
    M.opts.hyperparams = 1;  
    M.opts.hypertune   = 1; 
    M.opts.fsd         = 0;        
    M.opts.inner_loop  = 1;
    
    M.opts.objective   = 'gauss';%'qrmse_g';%'gauss';
    M.opts.criterion   = -inf;
    
    M.opts.factorise_gradients = 0;
    M.opts.normalise_gradients = 0;
    
    M.opts.memory_optimise = 0;
    M.opts.rungekutta      = 6;
    M.opts.bayesoptls      = 0;
    M.opts.updateQ         = 1; % do a grd ascent on Q but also weight by residual
    M.opts.crit            = [0 0 0 0];
    
    M.opts.userplotfun = @aodcmplotfun;
    
    M.opts.isNewton      = 0;
    M.opts.isQuasiNewton = 0;
    M.opts.isNewtonReg   = 0;      
    M.opts.isGaussNewton = 0;
    M.opts.isTrust       = 0;
    
    % order of dfdx: grads or curv & whether to orthogoanlise
    M.opts.order         = 1;
    M.opts.orthogradient = 1;
    M.opts.gradtol       = 1e-8;
        
    M.default_optimise([9],[18])
    
    M.update_parameters(M.Ep);

    M.default_optimise([9],[18])

    % save in DCM structures after optim 
    %----------------------------------------------------------------------
    DCM.M.pE = ppE;
    DCM.Ep = spm_unvec(M.Ep,DCM.M.pE);
    DCM.Cp = atcm.fun.reembedreducedcovariancematrix(DCM,M.CP);
    DCM.Cp = makeposdef(DCM.Cp);
    DCM.F = M.F;
    save(DCM.name); close all; clear global;
    
end

end
