function Setup(s)
% Top level script for setup and running TCM model
%
%
% AS

% cluster version: change dir
cd /cubric/scratch/sapas10/100B_VIS/

% Paths
addpath(genpath('~/spm12'));
addpath('/cubric/scratch/sapas10/tcm')

% Data & Design
%--------------------------------------------------------------------------
Data.Datasets     = 'VIS100_LIST.txt';  % textfile list of LFP SPM datasets (.txt)
Data.Design.X     = [];              % std/dev
Data.Design.name  = {'undefined'};         % condition names
Data.Design.tCode = [1];             % condition codes in SPM
Data.Design.Ic    = [1];             % channel indices
Data.Design.Sname = {'V1'};         % channel (node) names
Data.Prefix       = 'TCM_';      % outputted DCM prefix
Data.Datasets     = atcm.fun.ReadDatasets(Data.Datasets);

% Model space - T = ns x ns, where 1 = Fwd, 2 = Bkw
%--------------------------------------------------------------------------
T = [...
    0];
F = (T==1);
B = (T==2);
C = [1]';      % inputs
L = sparse(1,1); 

% Set up, over subjects
%--------------------------------------------------------------------------
%for s = 1:length(Data.Datasets)
    
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
    DCM.M.f  = @atcm.tc_dev;
    DCM.M.IS = @atcm.integrate;
    DCM.options.SpecFun = @atcm.fun.Afft;
    
    % Print Progress
    %----------------------------------------------------------------------
    fprintf('Running Dataset %d / %d\n',s,length(Data.Datasets));
    
    % Prepare Data
    %----------------------------------------------------------------------
    DCM.M.U            = sparse(diag(ones(Ns,1)));  %... ignore [modes]
    DCM.options.trials = tCode;                     %... trial code [GroupDataLocs]
    DCM.options.Tdcm   = [300 1000];                   %... peristimulus time
    DCM.options.Fdcm   = [4 80];                    %... frequency window
    DCM.options.D      = 1;                         %... downsample
    DCM.options.han    = 1;                         %... apply hanning window
    DCM.options.h      = 4;                         %... number of confounds (DCT)
    DCM.options.DoData = 1;                         %... leave on [custom]
    %DCM.options.baseTdcm   = [-200 0];                  %... baseline times [new!]
    DCM.options.Fltdcm = [4 80];                    %... bp filter [new!]

    DCM.options.analysis      = 'CSD';              %... analyse type
    DCM.xY.modality           = 'LFP';              %... ECD or LFP data? [LFP]
    DCM.options.spatial       = 'LFP';              %... spatial model [LFP]
    DCM.options.model         = 'tc6';              %... neural model
    DCM.options.Nmodes        = length(DCM.M.U);    %... number of modes

    % Alex additions - 1010 = use aFFT
    DCM.options.UseWelch      = 1010;
    DCM.options.FFTSmooth     = 2;
    DCM.options.UseButterband = [4 80];
    DCM.options.BeRobust      = 0;
    
    DCM.xY.name = DCM.Sname;
    DCM = atcm.fun.prepcsd(DCM);
    DCM.options.DATA = 1 ;      
    
    % Subfunctions
    %----------------------------------------------------------------------
    DCM = atcm.parameters(DCM,Ns);       % gets latet priors for tc nmm 
    DCM = atcm.complete(DCM);            % complete network specification
    
    % use the refined DEXPRO priors from now on:
    %load('AugSpectralPriors.mat');
    %DCM.M.pE = PX.pE;
    %DCM.M.pC = PX.pC;
    
    % New precision matrix
    %----------------------------------------------------------------------
    nf = length( DCM.options.Fdcm(1):DCM.options.Fdcm(end) );
    w  = DCM.options.Fdcm(1):DCM.options.Fdcm(end) ;
    y  = DCM.xY.y;
    n  = size(y{1},1);
    m  = size(y{1},2)*size(y{1},3);
    q  = spm_Q(1/2,n,1);
    q  = diag( .25+(w/w(end)) );   
    Q  = kron(speye(m,m),q);

    % peaks in data as high-precision points of interest
    %----------------------------------------------------------------------
    y0 = DCM.xY.y{1};
    [PKS,LOCS] = findpeaks(y0,w,'NPeaks',4);
    
    % make fitting the peaks twice as important as any other feature
    %----------------------------------------------------------------------
    for l = 1:length(LOCS)
        Q(LOCS(l),LOCS(l))   = Q(LOCS(l),LOCS(l)) + 2;
    end
    DCM.xY.Q  = (Q);
    DCM.xY.X0 = sparse(size(Q,1),0);

    % Fit the model using a non-Bayesian optimsation routine
    %DCM = atcm.optim.ainvert(DCM,'abc',length(find(spm_vec(DCM.M.pC))));
    
    DCM.M.nograph = 1;
    
    % DCM invert first
    [Qp,Cp,Eh,F] = atcm.optim.spm_nlsi_GN_as(DCM.M,DCM.xU,DCM.xY);
    DCM.M.pE = Qp;
    
    % now curvature inversion (AO.m)
    DCM = atcm.optim.ainvert(DCM,'ao');
    
    
    % Or, Fit the model using DCM
    %DCM = atcm.optim.dcminvert(DCM);
    %close; drawnow;
    
    %AllModels(s) = DCM;
    
    
%end
