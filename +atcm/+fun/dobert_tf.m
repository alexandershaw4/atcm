function DCM = dobert_tf(DCM)

% Set defaults and Get D filename
%-------------------------------------------------------------------------
try
    Dfile = DCM.xY.Dfile;
catch
    errordlg('Please specify data and trials');
    error('')
end

% ensure spatial modes have been computed (see spm_dcm_ssr)
%-------------------------------------------------------------------------
try
    DCM.M.U;
catch
    errordlg('Please estimate this model first');
    error('')
end

% load D
%--------------------------------------------------------------------------
try
    D = spm_eeg_load(Dfile);
catch
    try
        [p,f]        = fileparts(Dfile);
        D            = spm_eeg_load(f);
        DCM.xY.Dfile = fullfile(pwd,f);
    catch
        try
            [f,p]        = uigetfile('*.mat','please select data file');
            name         = fullfile(p,f);
            D            = spm_eeg_load(name);
            DCM.xY.Dfile = fullfile(name);
        catch
            warndlg([Dfile ' could not be found'])
            return
        end
    end
end

 
% indices of EEG channel (excluding bad channels)
%--------------------------------------------------------------------------
if ~isfield(DCM.xY, 'modality')
    [mod, list] = modality(D, 0, 1);

    if isequal(mod, 'Multimodal')
        qstr = 'Only one modality can be modelled at a time. Please select.';
        if numel(list) < 4
            % Nice looking dialog. Will usually be OK
            options = [];
            options.Default = list{1};
            options.Interpreter = 'none';
            DCM.xY.modality = questdlg(qstr, 'Select modality', list{:}, options);
        else
            % Ugly but can accomodate more buttons
            ind = menu(qstr, list);
            DCM.xY.modality = list{ind};
        end
    else
        DCM.xY.modality = mod;
    end
end



if ~isfield(DCM.xY, 'Ic')
    DCM.xY.Ic  = setdiff(D.meegchannels(DCM.xY.modality), D.badchannels);
else
    fprintf('Channel subspace selection: Using channel %d\n',DCM.xY.Ic);
end

Ic        = DCM.xY.Ic;
Nc        = length(Ic);
Nm        = size(DCM.M.U,2);
DCM.xY.Ic = Ic;

% options
%--------------------------------------------------------------------------
try
    DT    = DCM.options.D;
catch
    DT    = 1;
end
try
    trial = DCM.options.trials;
catch
    trial = D.nconditions;
end

% KRISH HAS COMMENTED OUT THIS CODE AS I WANT TO SPECIFY MY OWN DOWNSAMPLE RATE - THIS ACTUALLY HAS A BIG EFFECT ON THE SPECTRAL SHAPE
%
% check data are not oversampled (< 4ms)
%--------------------------------------------------------------------------
%if DT/D.fsample < 0.004
%    DT            = ceil(0.004*D.fsample);
%    DCM.options.D = DT;
%end


% get peristimulus times
%--------------------------------------------------------------------------
%try
    
    % time window and bins for modelling
    %----------------------------------------------------------------------
    DCM.xY.Time = 1000*D.time; % ms
    T1          = DCM.options.Tdcm(1);
    T2          = DCM.options.Tdcm(2);
    %[i, T1]     = min(abs(DCM.xY.Time - T1));
    %[i, T2]     = min(abs(DCM.xY.Time - T2));
    
    T1 = atcm.fun.findthenearest(T1, DCM.xY.Time);
    T2 = atcm.fun.findthenearest(T2, DCM.xY.Time);
    
    %KRISH HACK - SEE IF BASELINE HAS BEEN SPECIFIED
    DOBASE=0;
    try
           baseT1          = DCM.options.baseTdcm(1);
           baseT2          = DCM.options.baseTdcm(2);
           [i, baseT1]     = min(abs(DCM.xY.Time - baseT1));
           [i, baseT2]     = min(abs(DCM.xY.Time - baseT2));
           baseIt=[baseT1:DT:baseT2];
           DCM.xY.baseIt=baseIt;
           DOBASE=1;
    end
    
    % Time [ms] of down-sampled data
    %----------------------------------------------------------------------
    It          = [T1:DT:T2]';               % indices - bins
    DCM.xY.pst  = DCM.xY.Time(It);           % PST
    DCM.xY.It   = It;                        % Indices of time bins
    DCM.xY.dt   = DT/D.fsample;              % sampling in seconds
    Nb          = length(It);                % number of bins
    fprintf(1,'Time Indices=%d:%d in steps of %d\n',T1,T2,DT);
 
% get frequency range
%--------------------------------------------------------------------------
%try
    Hz1     = DCM.options.Fdcm(1);          % lower frequency
    Hz2     = DCM.options.Fdcm(2);          % upper frequency
    
    % Krish Hack to allow finer frequency resolution (mostly NOT for DCM
    % but for using this function to make CSDs...
    Fstep=1;
    try
        Fstep=DCM.options.FrequencyStep;
        fprintf(1,'Using fstep=%f\n',Fstep);
    catch
        Fstep=1;
    end

 
% Frequencies
%--------------------------------------------------------------------------
DCM.xY.Hz  = fix(Hz1:Fstep:Hz2);             % Frequencies
DCM.xY.Hz  = (Hz1:Fstep:Hz2);             % Frequencies: KRISH HACK TO ALLOW HIGHER RESOLUTION
Nf         = length(DCM.xY.Hz);        % number of frequencies
Ne         = length(trial);            % number of ERPs
 
 
% Cross spectral density for each trial type
%==========================================================================
condlabels = D.condlist;

try
    DCM.xY.y=DCM.Kinitial.ForceSpectra;
    fprintf('Using user-specified spectra...\n');
    DCM.xY.U    = DCM.M.U;
    DCM.xY.code = condlabels(trial);
    return;
end

try 
    bpf = DCM.options.UseButterband;
    fprintf('Bandpass filtering %d to %d Hz (%d Hz)\n',min(bpf),max(bpf),1/DCM.xY.dt);
end

for i = 1:Ne;
   clear Pfull
    
    % trial indices
    %----------------------------------------------------------------------
    %c = D.pickconditions(condlabels{trial(i)});
    c = D.indtrial(condlabels(trial(i)), 'GOOD');
    
    % use only the first 512 trial
    %----------------------------------------------------------------------
    %try c = c(1:512); end
    Nt    = length(c);
    
    % Get data
    %----------------------------------------------------------------------
    P     = zeros(Nf,Nm,Nm);
    if DOBASE==1,
       Pbase=zeros(Nf,Nm,Nm);
    end
    Fs = 1000/(DCM.xY.Time(2)-DCM.xY.Time(1));
    fprintf('\n');
        
    % loop trials
    for j = 1:Nt
        
        % print progress
        if j > 1; fprintf(repmat('\b',[1,length(str)])); end
        str = sprintf('Computing for trial %d/%d of condition %d',j,Nt,i);
        fprintf(str);
        if j == Nt; fprintf('\n'); end
        
        % Get Y for this condition
        Y    = full(double(D(Ic,It,c(j))'));%*DCM.M.U));
        
        trldat(j,:) = Y;
    end


    % load the virtual sensor and run a timefrequency analysis
    cfg.baseline = 'relchange';
    cfg.sampletimes = DCM.xY.pst;
    cfg.fsample = 1./DCM.xY.dt;
    cfg.filterorder = 4;
    FoI = linspace(Hz1,Hz2,DCM.options.nw);

    MatDat = trldat;

    tf{i} = atcm.fun.bert_singlechannel(MatDat,cfg,FoI,[-1 0]);
    
    agram = double(tf{i}.agram);
    
    agram = atcm.fun.HighResMeanFilt(agram,1,4);
    
    DCM.xY.y{i} = double(agram);
    DCM.xY.pst = tf{i}.sampletimes;
    
    DCM.xY.U    = DCM.M.U;
    DCM.xY.code = condlabels(trial);
    
end
    
    
    
    
    
    
    
    
    