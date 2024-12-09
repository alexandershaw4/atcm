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
Data.Datasets     = 'NewSZ.txt';%'MeanSZDatasets.txt';%'AllSZNoMerge.txt'; % textfile list of LFP SPM datasets (.txt)
Data.Design.X     = [];                % design matrix
Data.Design.name  = {'undefined'};     % condition names
Data.Design.tCode = [1];               % condition codes in SPM
Data.Design.Ic    = [1];               % channel indices
Data.Design.Sname = {'V1'};            % channel (node) names
Data.Prefix       = 'MTCM_';      % outputted DCM prefix
Data.Datasets     = atcm.fun.ReadDatasets(Data.Datasets);

%Data.Datasets     = {'NEW_MeanDataset.mat'};

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
    DCM.M.f  = @atcm.tc_hilge2;  % TCM2024;               % model function handle
    DCM.M.IS = @atcm.integrate_1;            % Alex integrator/transfer function
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
    DCM.options.Tdcm   = [300 1300];                   %... peristimulus time
    DCM.options.Fdcm   = fq;                    %... frequency window
    DCM.options.D      = 1;                         %... downsample
    DCM.options.han    = 0;                         %... apply hanning window
    DCM.options.h      = 1;                         %... number of confounds (DCT)
    DCM.options.DoData = 1;                         %... leave on [custom]
    %DCM.options.baseTdcm   = [-200 0];             %... baseline times [new!]
    DCM.options.Fltdcm = fq;                    %... bp filter [new!]
    DCM.options.UseButterband = fq;

    DCM.options.RegressionFFT = 0;

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

    DCM.xY.y{:} = abs(DCM.xY.y{:});

    DCM.xY.y{:}  = atcm.fun.agauss_smooth(DCM.xY.y{:},1);

    DCM.xY.y{:} = DCM.xY.y{:}';

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
    DCM.M.sim.dt  = 1./1000;
    DCM.M.sim.pst = 1000*((0:DCM.M.sim.dt:(1)-DCM.M.sim.dt)');
    DCM.M.burnin  = 0;
    
    % Input is d.c
    DCM.M.InputType = 0;

    % Use a 2-point RK method for integration
    DCM.M.intmethod = 0;%45;%1;%1;%45;

    %DCM.M.IntMethod = 'ode45_2';

    
    % No hamming on spectrum
    DCM.M.DoHamming = 0;

    load([p '/newpoints3.mat'],'pE','pC')

    pE = spm_unvec(spm_vec(pE)*0,pE);

    pC.ID = pC.ID * 0;
    pC.T  = pC.T *0;
    
    pE.J = pE.J-1000;    
    pE.J(1:8) = log([.6 .8 .4 .6 .4 .6 .4 .4]);
    %pE.J([7 8]) = log(.4);
    %pE.J([2 4]) = log(.8);
    %pC.ID = pC.ID + 1/8;
    pE.L = 2;
    pC.a = pC.a*0;

    pC.a(2) = 1/8;

    pC.C = 1/8;

    pC.J(1:8)=1/8;

    pC.d(1) = 1/8;
    
    pC.ID = ones(1,8)/8;
    %pC.S = ones(1,8)/8;

    %pE.f = [0 0];;
    %pC.f = [1 1]./8;

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

    

    % xx = load([p '/newx.mat']); DCM.M.x = spm_unvec(xx.x,DCM.M.x);
    % load('init_14dec','x');
    % DCM.M.x = spm_unvec(x,DCM.M.x);
    % 
    %x = atcm.fun.alexfixed(DCM.M.pE,DCM.M,1e-10);
    load('new_fp','x')
    DCM.M.x = spm_unvec(x,DCM.M.x);

    % precompute J and put in J which flags for the rhs to compute the
    % *delayed* update step
    [dx,J] = DCM.M.f(DCM.M.x,0,DCM.M.pE,DCM.M);
    DCM.M.J = J;

    % 
    % norm(DCM.M.f(DCM.M.x,0,DCM.M.pE,DCM.M))
    % 
    % % update Jacobian obtained from fixed point
    % [~,J] = DCM.M.f(DCM.M.x,0,DCM.M.pE,DCM.M);
    % DCM.M.J = J;
    % 
    %DCM.M.x = DCM.M.x*0;
    %DCM.M.x(:,:,1) = -70;

    fprintf('Finished...\n');

    % evaluate model and plot at fp
    %[y,X,pst] = simpleint(DCM.M.pE,DCM.M,DCM.xU);
    %subplot(211),plot(pst,X); subplot(212), plot(DCM.xY.Hz,y)
      
    fprintf('--------------- PARAM ESTIMATION (neural) ---------------\n');
    %fprintf('iteration %d\n',j);   

    
    % Construct an AO optimisation object
    M = AODCM(DCM);

    %M.bayesopt

    [parts,moments]=iterate_gauss(DCM.xY.y{:},2);
    for ii = 1:size(parts,1); QQ{ii} = diag(parts(ii,:)); end
    M.opts.Q = QQ;

    % Optimisation option set 1.
    M.opts.WeightByProbability = 0;

    M.opts.ismimo      = 1;     
    M.opts.doparallel  = 1;    

    M.opts.hyperparams = 1; 
    M.opts.hypertune   = 0; 
    M.opts.fsd         = 0;        
    M.opts.inner_loop  = 1;

    %M.opts.objective   = 'gaussfe';%_trace';%fe';%gauss_trace';%'gauss';%_trace';%'qrmse_g';%'gauss';
    M.opts.objective   = 'sse';

    M.opts.criterion   = -inf;

    M.opts.factorise_gradients = 0;
    M.opts.normalise_gradients = 1;

    M.opts.rungekutta      = 5;
    M.opts.surrls          = 0;
    M.opts.dopowell        = 0;
    M.opts.wolfelinesearch = 0;
    M.opts.bayesoptls      = 0; 
    M.opts.agproptls       = 0;
    M.opts.updateQ         = 0; 
    M.opts.crit            = [0 0 0 0];

    M.opts.userplotfun = @aodcmplotfun;

    % order of dfdx: grads or curv & whether to orthogoanlise
    M.opts.order         = 1;
    M.opts.gradtol       = 1e-8;

    M.default_optimise([1],[20]);
     
     M.update_parameters(M.Ep);

     %M.rungekutteopt(32)

     %M.update_parameters(M.Ep);
     
     M.default_optimise(1,20);

   

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
