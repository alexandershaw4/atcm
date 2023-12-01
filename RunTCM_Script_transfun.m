function RunTCM_Script_transfun(i)
% Top level script showing how to apply the thalamo-cortical neural mass
% model decribed in Shaw et al 2020 NeuroImage, to M/EEG data.
%
% This version using a linearisation and transfer function (numerical
% Laplace) rather than integration.
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
Data.Prefix       = 'stf_TCM_';      % outputted DCM prefix
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
    DCM.M.IS = @atcm.fun.alex_tf;            % Alex integrator/transfer function
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
        
    % Subfunctions and default priors
    %----------------------------------------------------------------------
    DCM = atcm.parameters(DCM,Ns);
    
    DCM.xY.y{:} = abs(DCM.xY.y{:});
    w = DCM.xY.Hz;

    %DCM.xY.y{:} = atcm.fun.awinsmooth(DCM.xY.y{:},2)';
        
    % other model options
    %----------------------------------------------------------------------
    DCM.M.solvefixed=0;      % oscillations == no fixed point search
    DCM.M.x = zeros(1,8,7);  % init state space: ns x np x nstates
    DCM.M.x(:,:,1)=-70;      % init pop membrane pot [mV]
        
    load('newpoints3','pE','pC')

    pC.d(1) = 1/8;
    pE.J(1:8) = log([.6 .8 .4 .6 .4 .6 .4 .4]);

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

    xx = load('newx'); DCM.M.x = spm_unvec(xx.x,DCM.M.x);
    x = atcm.fun.alexfixed(DCM.M.pE,DCM.M);
    DCM.M.x = spm_unvec(x,DCM.M.x);

    %DCM.M.x = FindSteadyState(DCM,20,1/1200); close;drawnow;
    %Y0 = spm_vec(feval(DCM.M.IS,DCM.M.pE,DCM.M,DCM.xU));
    %DCM.M.Y0 = Y0;
    fprintf('Finished...\n');
    
      
    fprintf('--------------- PARAM ESTIMATION ---------------\n');
    %fprintf('iteration %d\n',j);

    DCM.M.nograph = 1;

    [Qp,Cp,Eh,F] = spm_nlsi_GN(DCM.M,DCM.xU,DCM.xY);


    % % Construct an AO optimisation object
    % M = AODCM(DCM);
    % 
    % %QQ = spm_Q(1/2,length(w),1)*diag(w)*spm_Q(1/2,length(w),1)';
    % 
    % %M.opts.Q = QQ;
    % 
    % % Optimisation option set 1.
    % M.opts.EnforcePriorProb    = 0; 
    % M.opts.WeightByProbability = 0;
    % 
    % M.opts.ismimo      = 1;     
    % M.opts.doparallel  = 0;    
    % 
    % M.opts.hyperparams = 0; 
    % M.opts.ahyper      = 0;
    % M.opts.ahyper_p    = 0;
    % 
    % M.opts.hypertune   = 0; 
    % M.opts.fsd         = 0;        
    % M.opts.inner_loop  = 1;
    % 
    % M.opts.objective   = 'gauss_trace';%'gauss';%_trace';%'qrmse_g';%'gauss';
    % M.opts.criterion   = -inf;
    % 
    % M.opts.factorise_gradients = 0;
    % M.opts.normalise_gradients = 0;
    % 
    % M.opts.memory_optimise = 0;
    % M.opts.rungekutta      = 8;
    % M.opts.dopowell        = 0;
    % M.opts.wolfelinesearch = 0;
    % M.opts.bayesoptls      = 0;
    % M.opts.crit            = [0 0 0 0];
    % 
    % M.opts.isNewton      = 0;
    % M.opts.isQuasiNewton = 0;
    % M.opts.isNewtonReg   = 0;      
    % M.opts.isGaussNewton = 0;
    % M.opts.isTrust       = 0;
    % 
    % % order of dfdx: grads or curv & whether to orthogoanlise
    % M.opts.order         = 1;
    % M.opts.orthogradient = 1;
    % M.opts.gradtol       = 1e-8;
    % 
    % M.default_optimise([1],[20]);
    % 
    % M.update_parameters(M.Ep);
    % 
    % M.default_optimise(1,8);


    % M.update_parameters(M.Ep);
    %M.default_optimise([1],[8])

    % save in DCM structures after optim 
    %----------------------------------------------------------------------
    DCM.M.pE = ppE;
    DCM.Ep = Qp;%spm_unvec(M.Ep,DCM.M.pE);
    DCM.Cp = Cp;
    
    %DCM.Cp = atcm.fun.reembedreducedcovariancematrix(DCM,M.CP);
    %DCM.Cp = makeposdef(DCM.Cp);
    DCM.F = F;
    save(DCM.name); close all; clear global;
    
end

end
