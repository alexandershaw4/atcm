function NewRunVisGam_2023(i)
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

 
% Data & Design
%--------------------------------------------------------------------------
Data.Datasets     = 'MeanSZDatasets.txt';
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
    %---------------------------------------------------------------------
    DCM.M.f  = @atcm.experimental_models.spm_fx_tcm;               % model function handle
    DCM.M.IS = 'spm_csd_mtf';            % Alex integrator/transfer function
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
            
    % If using DCM inversion, select whether to block graph or not
    DCM.M.nograph = 0;
            
    % Feature function for the integrator [NOT USED]
    %----------------------------------------------------------------------
    DCM = atcm.complete(DCM);

    atcm.optim.dcminvert(DCM)
    
    
end