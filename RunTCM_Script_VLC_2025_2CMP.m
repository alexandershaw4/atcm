function RunTCM_Script_VLC_2025_2CMP(i)
% Multi-compartment version of TCM with soma + dendrite compoartns and
% separate presynaptic and postsynaptic 
%
% Top level script showing how to apply the thalamo-cortical neural mass
% model decribed in Shaw et al 2020 NeuroImage, to M/EEG data.
%
% This version using a linearisation and transfer function (numerical
% Laplace) rather than brute numerical integration.
%
% Requires atcm (thalamo cortical modelling package) and aoptim
% (optimisation package)
%
% atcm: https://github.com/alexandershaw4/atcm
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
Data.Design.Sname = {'PBVE'};            % channel (node) names
Data.Prefix       = 'nVLC_TFD_TCM_';      % outputted DCM prefix
Data.Datasets     = atcm.fun.ReadDatasets(Data.Datasets);

% Model space - T = ns x ns, where 1 = Fwd, 2 = Bkw
%--------------------------------------------------------------------------
T = [... % this is a 1-node model; nothing to put here...
    0];
F = (T==1);
B = (T==2);
C = [1]';          % input(s)
L = sparse(1,1);

[p]=fileparts(which('atcm.integrate_1'));p=strrep(p,'+atcm','');addpath(p);


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
    %DCM.M.f  = @atcm.tc_hilge2;               % model function handle
    DCM.M.f = @atcm.tc_twocmp_stp;
    DCM.M.IS = @atcm.fun.Alex_LaplaceTFwD;            % Alex integrator/transfer function
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
    DCM.options.Tdcm   = [1 2000];                   %... peristimulus time
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

    DCM.xY.y{:}  = agauss_smooth(abs(DCM.xY.y{:}),1)';
        
    % Subfunctions and default priors
    %----------------------------------------------------------------------
    DCM = atcm.parameters(DCM,Ns);
            
    % other model options
    %----------------------------------------------------------------------
    DCM.M.solvefixed=0;      % 
    DCM.M.x = zeros(1,8,7);  % init state space: ns x np x nstates
    DCM.M.x(:,:,1)=-70;      % init pop membrane pot [mV]

    % --- Expand state vector: 1×8×7  ->  1×8×10
    % Indices (for clarity)
    iVs = 1; iGE = 2; iGI = 3; iGN = 4; iGB = 5; iGM = 6; iGH = 7;
    iVd = 8; iR  = 9; iU  = 10;

    % If M.x existed with 7 states, expand to 10 (Vd, R, uSTP)
    if size(DCM.M.x,3) < 10
        Xnew = zeros(size(DCM.M.x,1), size(DCM.M.x,2), 10);
        Xnew(:,:,1:size(DCM.M.x,3)) = DCM.M.x;

        % Initialise new states
        % Vd starts near Vs; R~0.9 resources available; uSTP ~ baseline use
        Xnew(:,:,iVd) = Xnew(:,:,iVs);       % dendrite voltage ~ soma
        Xnew(:,:,iR)  = 0.9;                 % resources near 1
        Xnew(:,:,iU)  = 0.2;                 % 'use' (will be overridden by U0 prior below)

        DCM.M.x = Xnew;
    end

        
    load([p '/newpoints3.mat'],'pE','pC')

    pE = spm_unvec(spm_vec(pE)*0,pE);

    pC.ID = pC.ID * 0;
    pC.T  = pC.T *0;
    
    pE.J = DCM.M.x(:)*0-1000;    
    %pE.J(1:8) = log([.6 .8 .4 .6 .4 .6 .4 .4]);
    pE.J(1:8) = log([.2 .99 .1 .8 .1 .2 .05 .1]);
    pC.J = pE.J*0;
    %pC.ID = pC.ID + 1/8;
    pE.L = 0;
    pC.a = pC.a*0;

    pE.Gb = pE.H;
    pC.Gb = [1   0   0   0   0   0   0   0;
             0   1   1   0   0   0   0   0;
             0   0   1   0   0   0   0   0;
             0   0   0   1   1   0   0   0;
             0   0   0   0   1   0   0   0;
             0   0   0   0   1   1   0   0;
             0   0   0   0   0   0   0   0;
             0   0   0   0   0   0   1   0]/64;

    %pC.J(1:8)=1/8;
    pC.d(1) = 1/8;
    pC.d(3) = 1/8;

    logit = @(p) log(p./(1-p));

    % === STP (per presyn population: 8×1 each, used across its outputs)
    if ~isfield(pE,'prel'),  pE.prel  = log(0.6)*ones(8,1); end   % baseline release prob (log-space)
    if ~isfield(pE,'tauR'),  pE.tauR  = 0.6*ones(8,1);      end   % R recovery time (s)
    if ~isfield(pE,'tauU'),  pE.tauU  = 0.2*ones(8,1);      end   % u facilitation time (s)
    if ~isfield(pE,'U0'),    pE.U0    = logit(0.2)*ones(8,1); end % baseline 'use' (stored in logit; code does logistic)

    if ~isfield(pC,'prel'),  pC.prel  = ones(8,1)/8; end   % baseline release prob (log-space)
    if ~isfield(pC,'tauR'),  pC.tauR  = ones(8,1)/8;      end   % R recovery time (s)
    if ~isfield(pC,'tauU'),  pC.tauU  = ones(8,1)/8;      end   % u facilitation time (s)
    if ~isfield(pC,'U0'),    pC.U0    = ones(8,1)/8; end % baseline 'use' (stored in logit; code does logistic)


    % === Two-compartment routing & coupling
    if ~isfield(pE,'gc'),      pE.gc      = log(3);       end     % axial coupling (log-space; ~3 mS default)
    if ~isfield(pE,'w_dend'),  pE.w_dend  = [0.8 0.9];   end     % [AMPA_to_dend, NMDA_to_dend] in [0,1]
    if ~isfield(pE,'w_soma'),  pE.w_soma  = [0.9 0.9];   end     % [GABAa_to_soma, GABAb_to_soma] in [0,1]
    
    if ~isfield(pC,'gc'),      pC.gc      = 1/8;       end     % axial coupling (log-space; ~3 mS default)
    if ~isfield(pC,'w_dend'),  pC.w_dend  = [1 1]/8;   end     % [AMPA_to_dend, NMDA_to_dend] in [0,1]
    if ~isfield(pC,'w_soma'),  pC.w_soma  = [1 1]/2;   end     % [GABAa_to_soma, GABAb_to_soma] in [0,1]

    pE.scale = zeros(4,1);
    pC.scale = ones(4,1)/8;

    % PC0.prel   = (0.3^2) * ones(8,1);   % log-space SD ~0.3
    % PC0.tauR   = (0.2^2) * ones(8,1);
    % PC0.tauU   = (0.1^2) * ones(8,1);
    % PC0.U0     = (0.5^2) * ones(8,1);   % logit-space SD
    % PC0.gc     = (0.3^2);               % log-space SD
    % PC0.w_dend = (0.2^2) * ones(1,2);   % direct [0,1] (clamped in model)
    % PC0.w_soma = (0.2^2) * ones(1,2);


    % Make changes here;
    %-----------------------------------------------------------
   
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

    DCM.M.x(:,:,iU) = repmat( 1./(1+exp(-DCM.M.pE.U0(:)))', size(DCM.M.x,1), 1 );  % logistic(U0)

    %xx = load([p '/newx.mat']); 
    %DCM.M.x = spm_unvec(xx.x,DCM.M.x);
    %load('init_14dec','x');
    %DCM.M.x = spm_unvec(x,DCM.M.x);

    x = atcm.fun.alexfixed(DCM.M.pE,DCM.M,1e-10);
    DCM.M.x = spm_unvec(x,DCM.M.x);

    norm(DCM.M.f(DCM.M.x,0,DCM.M.pE,DCM.M))

    fprintf('Finished...\n');
    
          
    fprintf('--------------- PARAM ESTIMATION ---------------\n');
    %fprintf('iteration %d\n',j);

    % Fit with DCM VB routine:
    %----------------------------------------------------------------------
    %[Qp,Cp,Eh,F] = spm_nlsi_GN(DCM.M,DCM.xU,DCM.xY);
    

    % Fit with LM (Log Likelihood estimation):
    %----------------------------------------------------------------------
    M = aFitDCM(DCM)

    M.aloglikVLtherm([],0);
    M.update_parameters(M.Ep)

    [y,w,G,s] = feval(DCM.M.IS,spm_unvec(M.Ep,DCM.M.pE),DCM.M,DCM.xU);

    numit = 0;
    while cdist(DCM.xY.y{:}',y{:}') > (1/2) && numit < 8
        numit = numit + 1;
        %M.aloglik;
        M.aloglikVLtherm([],0);
        %M.aloglikFE;
        M.update_parameters(M.Ep);
        [y,w,G,s] = feval(DCM.M.IS,spm_unvec(M.Ep,DCM.M.pE),DCM.M,DCM.xU);
    end
    
    Qp = spm_unvec(M.Ep,DCM.M.pE);
    Cp = M.CP;
    F  = M.F;


    % Fit with LM (Free energy estimation):
    %----------------------------------------------------------------------
    % M = aFitDCM(DCM)
    % 
    % M.aloglikFE
    % M.update_parameters(M.Ep)
    % M.aloglikFE
    % M.update_parameters(M.Ep)
    % M.aloglikFE
    %Qp = spm_unvec(M.Ep,DCM.M.pE);
    %Cp = M.CP;


    % save in DCM structures after optim 
    %----------------------------------------------------------------------
    DCM.M.pE = ppE;
    DCM.Ep = Qp;%spm_unvec(M.Ep,DCM.M.pE);
    DCM.Cp = Cp;

    DCM.M.sim.dt  = 1./600;
    DCM.M.sim.pst = 1000*((0:DCM.M.sim.dt:(2)-DCM.M.sim.dt)');

    [y,w,G,s] = feval(DCM.M.IS,DCM.Ep,DCM.M,DCM.xU);

    DCM.pred = y;
    DCM.w = w;
    DCM.G = G;
    DCM.series = s;
    
    %DCM.Cp = atcm.fun.reembedreducedcovariancematrix(DCM,M.CP);
    %DCM.Cp = makeposdef(DCM.Cp);
    DCM.F  = F;%M.FreeEnergyF;
    %DCM.Cp = M.CP;
    save(DCM.name); close all; clear global;
    
end

end
