function NewRunDEPKET_2022(i)
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
clear global;

% for cubric system, add paths -
% addpath(genpath('~/spm12'));
% addpath(genpath('/home/sapas10/code/atcm/'));
% addpath(genpath('/home/sapas10/code/aoptim/'));
 
% Data & Design
%--------------------------------------------------------------------------
Data.Datasets     = 'ketdep_list.txt';
Data.Design.X     = [];                % design matrix
Data.Design.name  = {'undefined'};     % condition names
Data.Design.tCode = [1];               % condition codes in SPM
Data.Design.Ic    = [1];               % channel indices
Data.Design.Sname = {'V1'};            % channel (node) names
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
    
    % Frequency range of interest
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

    % 1010 == use atcm.fun.AFFT.m
    DCM.options.UseWelch      = 1010;
    DCM.options.FFTSmooth     = 3;
    %DCM.options.UseButterband = fq;
    DCM.options.BeRobust      = 0;
    DCM.options.FrequencyStep = 1;     
            
    DCM.xY.name = DCM.Sname;
    DCM = atcm.fun.prepcsd(DCM);
    DCM.options.DATA = 1 ;      
                
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
    Nf = length(w);
    Q  = {spm_Q(1/2,Nf,1)*diag(DCM.M.Hz)*spm_Q(1/2,Nf,1)};
    Qw = Qw * Q{:};

    % Newton-Cotes integration parameters
    %----------------------------------------------------------------------
    DCM.M.sim.dt  = 1./300;
    DCM.M.sim.pst = 1000*((0:DCM.M.sim.dt:(3)-DCM.M.sim.dt)');
    DCM.M.burnin  = 300;
    DCM.M.intmethod = 44;
    
    % Input is an ERP
    DCM.M.InputType = 0;
    DCM.M.pE.C = log(.01);
    
    % only interested in real psd rn
    %----------------------------------------------------------------------
    DCM.xY.y{1} = real(DCM.xY.y{1});
    DCM.M.y     = DCM.xY.y;
    
    ppE = DCM.M.pE;
    
    % Optimise using AO.m                                                     
    %----------------------------------------------------------------------
    M = AODCM(DCM);

    % Bias and feature selection
    M.opts.Q  = real(Qw);  
    M.opts.FS = @(x) [real(sqrt(x))];
    M.opts.FS = @(x) [real( spm_vec(atcm.fun.Pf2VMD(x,3)) )];    
    
    % opt set 1.
    M.opts.EnforcePriorProb=0; % forcibly constrain parameters to within prior dist
    M.opts.ismimo=0;        % compute dfdp elementwise on vector-output function
    M.opts.doparallel=1;    % use parfor loops when poss, incl for df/dx
    M.opts.hyperparams=1;   % hyperparameter tuning
    M.opts.fsd = 1;         % fixed-step for derivatives
    M.opts.corrweight = 0;  % weight log evidence by correlation 
    M.opts.inner_loop = 3;
        
    M.opts.objective = 'qrmse'; % objective (error) function
    M.opts.criterion = 1e-3;
    
    %M.opts.isGaussNewton=1;
    
    M.default_optimise([7],[28])
    
    % afterward, use AODCM object to loop through the optimisation steps
    % for a visualisation:
    %-------------------------------------------------------------------
%     for i = 1:28; 
%        dydp(i,:) = spm_vec(M.opts.fun(M.history.p{i}));
%     end
    
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
    M.opts.hyperparams=1;   % hyperparameter tuning
    M.opts.fsd=0;           % fixed-step for derivatives
    M.opts.corrweight = 0;  % weight log evidence by correlation 
    
    M.opts.objective = 'qrmse'; % objective (error) function
    M.opts.criterion = 1e-3;

    M.default_optimise([7],[8])
      
    % reinstate the actual priors before saving
    DCM.M.pE = ppE;
    
    save(DCM.name); close; clear global;    
    
end