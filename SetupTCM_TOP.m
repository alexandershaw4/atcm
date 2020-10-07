% Top level script showing how to apply the thalamo-cortical neural mass
% model decribed in Shaw et al 2020 NeuroImage, to M/EEG data. 
%
% Besides the SPM built-ins, all the dependent functions are contained in
% the +atcm/ package.
%
% I include here examples of how to set up both single-region and
% multi-region models.
%
% As of late 2020, I am having more success using an alternative
% optimisation routine 'aoptim', available at https://github.com/alexandershaw4/aoptim
% This routine is similar to the DCM one (spm_nlsi)GN.m) - it still optimises 
% free energy using a gradient descent, but has some extra options such as 
% momentum parameters and line search. 
%
%
% AS2020

% EXAMPLE ONE NODE SETUP:
%==========================================================================

% Data & Design
%--------------------------------------------------------------------------
Data.Datasets     = 'KetDatasets.txt';  % textfile list of LFP SPM datasets (.txt)
Data.Design.X     = [];              % std/dev
Data.Design.name  = {'undefined'};         % condition names
Data.Design.tCode = [1];             % condition codes in SPM
Data.Design.Ic    = [1];             % channel indices
Data.Design.Sname = {'V1'};         % channel (node) names
Data.Prefix       = 'TCM_';      % outputted DCM prefix
Data.Datasets     = ReadverifyDatasets(Data.Datasets);

% Model space - T = ns x ns, where 1 = Fwd, 2 = Bkw
%--------------------------------------------------------------------------
T = [... % this is a 1-node model; nothing to put here...
    0];
F = (T==1);
B = (T==2);
C = [1]';          % input(s)
L = sparse(1,1); 


% EXAMPLE EIGHT NODE SETUP:
%==========================================================================

% % Data & Design
% %--------------------------------------------------------------------------
% Data.Datasets     = {'314_27934_VisAAL'};  % textfile list of LFP SPM datasets (.txt)
% Data.Design.X     = [];              % std/dev
% Data.Design.name  = {'Stim'};         % condition names
% Data.Design.tCode = [1];             % condition codes in SPM
% % l / r thal: 77/78
% % l / r calc: 43/44
% % l / r sup occ: 49/50 (dorsal stream)
% % l / r inf occ: 53/54 (central stream)
% Data.Design.Ic    = [77 78 43 44 49 50 53 54];          % channel indices
% Data.Design.Sname = {'LTh' 'RTh' 'LV1' 'RV1' 'LSupOc' 'RSupOv' 'LInfOc' 'RInfOc'};  % channel (node) names
% Data.Prefix       = 'TCM_';      % outputted DCM prefix
% 
% % Model space - T = ns x ns, where 1 = Fwd, 2 = Bkw
% %--------------------------------------------------------------------------
% T = [...
%     0 0 2 0 0 0 0 0;  % th
%     0 0 0 2 0 0 0 0;  % th
%     1 0 0 0 2 0 2 0;  % v1
%     0 1 0 0 0 2 0 2;  % v1
%     0 0 1 0 0 0 0 0;  % supoc
%     0 0 0 1 0 0 0 0;  % supoc 
%     0 0 1 0 0 0 0 0;  % infoc
%     0 0 0 1 0 0 0 0 ];% infoc
% 
% F = (T==1);
% B = (T==2);
% C = [1 1 0 0 0 0 0 0]';      % inputs
% L = sparse(8,8); 


% THE REST IS COMMON ACROSS SETUPS: 
%                                   DONT EDIT, APART FROM 'PREPARE DATA'
%--------------------------------------------------------------------------





% Set up, over subjects
%--------------------------------------------------------------------------
for s = 1:length(Data.Datasets)
    
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
    DCM.M.f  = @atcm.tc_dev_dev;               % model function handle
    DCM.M.IS = @atcm.integrate3;              % Alex integrator/transfer function
    %DCM.M.IS = 'spm_csd_mtf';              % default DCM transfer function
    DCM.options.SpecFun = @atcm.fun.Afft;  % fft function for IS
    
    % Print Progress
    %----------------------------------------------------------------------
    fprintf('Running Dataset %d / %d\n',s,length(Data.Datasets));
    
    % Prepare Data
    %----------------------------------------------------------------------
    DCM.M.U            = sparse(diag(ones(Ns,1)));  %... ignore [modes]
    DCM.options.trials = tCode;                     %... trial code [GroupDataLocs]
    DCM.options.Tdcm   = [300 1000];                   %... peristimulus time
    DCM.options.Fdcm   = [2 120];                    %... frequency window
    DCM.options.D      = 1;                         %... downsample
    DCM.options.han    = 1;                         %... apply hanning window
    DCM.options.h      = 4;                         %... number of confounds (DCT)
    DCM.options.DoData = 1;                         %... leave on [custom]
    %DCM.options.baseTdcm   = [-200 0];             %... baseline times [new!]
    DCM.options.Fltdcm = [2 120];                    %... bp filter [new!]

    DCM.options.analysis      = 'CSD';              %... analyse type
    DCM.xY.modality           = 'LFP';              %... ECD or LFP data? [LFP]
    DCM.options.spatial       = 'LFP';              %... spatial model [LFP]
    DCM.options.model         = 'tc6';              %... neural model
    DCM.options.Nmodes        = length(DCM.M.U);    %... number of modes

    % Alex additions - 1010 = use atcm.fun.AFFT.m
    DCM.options.UseWelch      = 1010;
    DCM.options.FFTSmooth     = 0;
    DCM.options.UseButterband = [2 120];
    DCM.options.BeRobust      = 0;
    DCM.options.FrequencyStep = 0.5;        % use .5 Hz steps
    
    DCM.xY.name = DCM.Sname;
    DCM = atcm.fun.prepcsd(DCM);
    DCM.options.DATA = 1 ;      
    
    % Subfunctions
    %----------------------------------------------------------------------
    DCM = atcm.parameters(DCM,Ns);       % gets latet priors for tc nmm     
    DCM = atcm.complete(DCM);            % complete network specification
            
    % if using AOPTIM for inversion, invoke the linear model g(x) output by
    % placing data (DCM.xY.y) in model struct - DCM.M.y
    DCM.M.y  = DCM.xY.y;
    DCM.M.Hz = DCM.M.Hz;
    
    % If using DCM inversion, select whether to block graph or not
    DCM.M.nograph = 0;
    
    % Final options for integrator
    DCM.M.fmethod = 'none';
    DCM.M.DoEnv   = 0;
    
    % Fit the model using DCM inversion routine:
    %-----------------------------------------------------
    %DCM = atcm.optim.dcminvert(DCM);
    %close; drawnow;
    
    % Or, invert (optimise) using aoptim (AO.m) - works better!
    %-----------------------------------------------------
    M = AODCM(DCM);
    M.default_optimise();
    
    % Extract the things we need from the optimisation object
    EP = M.Ep;
    F  = M.F;
    CP = M.CP;
    History = M.history;

    EP = spm_unvec( spm_vec(EP), DCM.M.pE);
    
    % re-embed the reduced covariance matrix into full-model space
    CP1 = atcm.fun.reembedreducedcovariancematrix(DCM,CP1);
    
    % save outputs when using AO.m
    save(DCM.name,'DCM','EP','F','CP','History','M');

    
    % n.b. 
    % EP = posteriors, F = (sign flipped) F-value, CP = coviarance of
    % (active_ parameters), History = history
    % structure from the optimisation routine
    
    
    
end
