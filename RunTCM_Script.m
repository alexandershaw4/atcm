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
Data.Datasets     = 'NewMeansSZ.txt';%'MeanSZDatasets.txt';%'AllSZNoMerge.txt'; % textfile list of LFP SPM datasets (.txt)
Data.Design.X     = [];                % design matrix
Data.Design.name  = {'undefined'};     % condition names
Data.Design.tCode = [1];               % condition codes in SPM
Data.Design.Ic    = [1];               % channel indices
Data.Design.Sname = {'V1'};            % channel (node) names
Data.Prefix       = 'Tscale_FP_TCM_';      % outputted DCM prefix
Data.Datasets     = atcm.fun.ReadDatasets(Data.Datasets);

[p]=fileparts(which('atcm.integrate_1'));p=strrep(p,'+atcm','');addpath(p);


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
    fq =  [3 90];
    
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
    DCM.options.FrequencyStep = 1;
    
    DCM.xY.name = DCM.Sname;
    DCM = atcm.fun.prepcsd(DCM);
    DCM.options.DATA = 1 ;


    DCM.xY.y{:}  = atcm.fun.agauss_smooth(DCM.xY.y{:},.6);

    %DCM.xY.y{:} = atcm.fun.awinsmooth(DCM.xY.y{:},2)';

    % also without the robust fitting to get the residual
    %DCMo = DCM;
    %DCMo.options.BeRobust=0;
    %DCMo = atcm.fun.prepcsd(DCMo);
    %r = DCMo.xY.y{1} - DCM.xY.y{1};

    % amount of smoothing scales linearly with frequency step
    %SmoothingK = 4./DCM.options.FrequencyStep;
        
    % Subfunctions and default priors
    %----------------------------------------------------------------------
    DCM = atcm.parameters(DCM,Ns);
    
    DCM.xY.y{:} = abs(DCM.xY.y{:});
    w = DCM.xY.Hz;

    %DCM.xY.y{:} = atcm.fun.awinsmooth(DCM.xY.y{:},2)';
    
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
    DCM.M.solvefixed=0;      % 
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
    DCM.M.intmethod = 45;%45;

    % No hamming on spectrum
    DCM.M.DoHamming = 0;

    load([p '/newpoints3.mat'],'pE','pC')

    pE = spm_unvec(spm_vec(pE)*0,pE);

    pC.ID = pC.ID * 0;
    pC.T  = pC.T *0;

    
    pE.J = pE.J-1000;
    pE.J([1 2 4]) = log([.4 .8 .6]);
    
    pE.J(1:8) = log([.6 .8 .4 .6 .4 .6 .4 .4]);
        
    pC.ID = pC.ID + 1/8;
    %pC.J(1:8) = 1/8;


    pE.L = 0;

    DCM.M.pE = pE;
    DCM.M.pC = pC;
    
    %DCM.M.dmd=1;

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
    %x = atcm.fun.alexfixed(DCM.M.pE,DCM.M,1e-10);
    load('init_14dec','x');
    DCM.M.x = spm_unvec(x,DCM.M.x);

    %x = atcm.fun.alexfixed(DCM.M.pE,DCM.M,1e-10);
    %DCM.M.x = spm_unvec(x,DCM.M.x);

    %xx0 = DCM.M.f(DCM.M.x,0,DCM.M.pE,DCM.M);
    %DCM.M.x = spm_unvec(xx0,DCM.M.x);

    %Y0 = spm_vec(feval(DCM.M.IS,DCM.M.pE,DCM.M,DCM.xU));
    %DCM.M.Y0 = Y0;
    
    fprintf('Finished...\n');

    %load('GausFFTMat','Mt');
    %DCM.M.GFFTM = Mt;
      
    fprintf('--------------- PARAM ESTIMATION (neural) ---------------\n');
    %fprintf('iteration %d\n',j);   

    % Construct an AO optimisation object
    M = AODCM(DCM);

    % Optimisation option set 1.
    M.opts.EnforcePriorProb    = 0; 
    M.opts.WeightByProbability = 0;

    M.opts.ismimo      = 1;     
    M.opts.doparallel  = 1;    

    M.opts.hyperparams = 1; 
    M.opts.ahyper      = 0;
    M.opts.ahyper_p    = 0;

    M.opts.hypertune   = 0; 
    M.opts.fsd         = 0;        
    M.opts.inner_loop  = 1;

    M.opts.objective   = 'gaussfe';%_trace';%fe';%gauss_trace';%'gauss';%_trace';%'qrmse_g';%'gauss';
    M.opts.criterion   = -inf;

    M.opts.factorise_gradients = 0;
    M.opts.normalise_gradients = 0;

    M.opts.memory_optimise = 0;
    M.opts.rungekutta      = 4;
    M.opts.surrls          = 0;
    M.opts.dopowell        = 0;
    M.opts.wolfelinesearch = 0;
    M.opts.bayesoptls      = 0;
    M.opts.agproptls       = 0;
    M.opts.updateQ         = 0; 
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

    M.default_optimise([1],[20]);
     
     M.update_parameters(M.Ep);
     
     M.default_optimise(1,20);

    % opts = AO('options');
    % 
    % % model and parameters and data
    % opts.fun = @(P,M) spm_vec(DCM.M.IS(spm_unvec(P,DCM.M.pE),DCM.M,DCM.xU));
    % 
    % opts.x0 = spm_vec(DCM.M.pE);
    % opts.M  = DCM.M;
    % opts.V  = spm_vec(DCM.M.pC);
    % opts.y  = spm_vec(DCM.xY.y);
    % 
    % %opts.isDynamicalSS = 1;
    % 
    % opts.maxit = 24;
    % 
    % opts.EnforcePriorProb    = 0; 
    % opts.WeightByProbability = 0;
    % 
    % opts.ismimo      = 1;     
    % opts.doparallel  = 1;    
    % 
    % opts.hyperparams = 1; 
    % opts.ahyper      = 0;
    % opts.ahyper_p    = 0;
    % 
    % opts.hypertune   = 0; 
    % opts.fsd         = 0;        
    % opts.inner_loop  = 1;
    % 
    % opts.objective   = 'sse';%_trace';%'gauss';%_trace';%'qrmse_g';%'gauss';
    % opts.criterion   = -inf;
    % 
    % opts.factorise_gradients = 0;
    % opts.normalise_gradients = 1;
    % 
    % opts.memory_optimise = 0;
    % opts.rungekutta      = 4;
    % opts.dopowell        = 0;
    % opts.wolfelinesearch = 0;
    % opts.bayesoptls      = 0;
    % opts.agproptls       = 0;
    % opts.updateQ         = 0; 
    % opts.crit            = [0 0 0 0];
    % 
    % opts.allow_worsen = 1000;
    % 
    % opts.isNewton      = 0;
    % opts.isQuasiNewton = 0;
    % opts.isNewtonReg   = 0;      
    % opts.isGaussNewton = 1;
    % opts.isTrust       = 0;
    % 
    % % order of dfdx: grads or curv & whether to orthogoanlise
    % opts.im            = 0;
    % opts.order         = 1;
    % opts.orthogradient = 1;
    % opts.gradtol       = 1e-8;
    % 
    % [X,F,CV,~,Hi] = AO(opts);    
    % % 

%     fprintf('--------------- STATE ESTIMATION ---------------\n');
%     fprintf('iteration %d\n',j);
% 
%     % one-iteration states search
%     opts     = AO('options');
%     opts.fun = @(x) spm_vec(feval(DCM.M.IS,DCM.M.pE,DCM.M,DCM.xU,x));
%     opts.x0  = DCM.M.x(:);
%     opts.V   = ones(56,1)/8;
%     opts.y   = DCM.xY.y{:};
% 
%     opts.EnforcePriorProb    = 0; 
%     opts.WeightByProbability = 0;
% 
%     opts.ismimo      = 1;     
%     opts.doparallel  = 1;    
% 
%     opts.hyperparams = 1; 
%     opts.ahyper      = 1;
%     opts.ahyper_p    = 1;
% 
%     opts.hypertune   = 1; 
%     opts.fsd         = 0;        
%     opts.inner_loop  = 1;
% 
%     opts.objective   = 'gauss_trace';%'gauss';%_trace';%'qrmse_g';%'gauss';
%     opts.criterion   = -inf;
% 
%     opts.factorise_gradients = 0;
%     opts.normalise_gradients = 0;
%     opts.memory_optimise = 0;
%     opts.rungekutta      = 8;
%     opts.dopowell        = 0;
%     opts.wolfelinesearch = 0;
%     opts.bayesoptls      = 0;
%     opts.updateQ         = 1; 
%     opts.crit            = [0 0 0 0];
% 
%     opts.isNewton      = 0;
%     opts.isQuasiNewton = 0;
%     opts.isNewtonReg   = 0;      
%     opts.isGaussNewton = 0;
%     opts.isTrust       = 0;
% 
%     % order of dfdx: grads or curv & whether to orthogoanlise
%     opts.order         = 1;
%     opts.orthogradient = 1;
%     opts.gradtol       = 1e-8;
% 
%     opts.maxit = 2;
% 
%     [X,F,Cp] = AO(opts);
% 
%     %update initial states ready for parameter estimation
%     DCM.M.x = spm_unvec(X,DCM.M.x);
% 
% end
   

    % M.update_parameters(M.Ep);
    %M.default_optimise([1],[8])

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
