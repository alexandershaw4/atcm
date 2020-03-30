function DCM = preptf(DCM)
% gets cross-spectral density data-features using a VAR model
% FORMAT DCM = spm_dcm_csd_data_as(DCM)
% DCM    -  DCM structure
% requires
%
%    DCM.xY.Dfile        - name of data file
%    DCM.M.U             - channel subspace
%    DCM.options.trials  - trial to evaluate
%    DCM.options.Tdcm    - time limits
%    DCM.options.Fdcm    - frequency limits
%    DCM.options.D       - Down-sampling
%
% sets
%
%    DCM.xY.pst     - Peristimulus Time [ms] sampled
%    DCM.xY.dt      - sampling in seconds [s] (down-sampled)
%    DCM.xY.U       - channel subspace
%    DCM.xY.y       - cross spectral density over sources
%    DCM.xY.csd     - cross spectral density over sources
%    DCM.xY.It      - Indices of time bins
%    DCM.xY.Ic      - Indices of good channels
%    DCM.xY.Hz      - Frequency bins
%    DCM.xY.code    - trial codes evaluated



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
        
    for j = 1:Nt
        
        % print progress
        if j > 1; fprintf(repmat('\b',[1,length(str)])); end
        str = sprintf('Computing for trial %d/%d of condition %d',j,Nt,i);
        fprintf(str);
        if j == Nt; fprintf('\n'); end
        
        % Get Y for this condition
        Y    = full(double(D(Ic,It,c(j))'));%*DCM.M.U));
        Ymod = Y;
        
        
        Fq = DCM.options.UseButterband;
        OLP = 16;
        
        % Convert edges to stepwise vector
        if length(Fq) == 2
            Fq = linspace(Fq(1),Fq(2), ( Fq(2)-Fq(1) )*4 );
        end
        
        clear fPf iPf
        
        % FORWARD FILTERS
        for fstep = 1:length(Fq)-OLP
            
            [dY,ix] = atcm.fun.padtimeseries(Y);
            fY = atcm.fun.bandpassfilter(dY,Fs,[Fq(fstep) Fq(fstep+OLP)]);
            fY = abs(hilbert(fY));
            fPf(fstep,j,:) = fY(ix);
        
        end  
        fPf(end+1:end+OLP,j,:) = repmat( fPf(end,j,:) ,[OLP 1] );
        
        % BACKWARD FILTERS
        iFq = fliplr(Fq);
        for fstep = 1:length(Fq)-OLP
            [dY,ix] = atcm.fun.padtimeseries(Y);
            fY = atcm.fun.bandpassfilter(dY,Fs,[iFq(fstep) iFq(fstep+OLP)]);
            fY = abs(hilbert(fY));
            iPf(fstep,j,:) = fY(ix);
        end     
        iPf(end+1:end+OLP,j,:) = repmat( iPf(end,j,:) ,[OLP 1] );
        
        % Flip backward matrix
        iPf = flipud(iPf);
        
        Pf = (fPf + iPf) / 2;
    end
    
    Pfull = (squeeze(mean(Pf,2))'./Fq.^2)' ;
        
        P = Pfull;
        
        P(find(isinf(P)))=0;
        P(find(isnan(P)))=0;
        %P(find(P<=0))=0;
        
        %P=P-min(P);
        
%         if isfield(DCM.options,'han') && DCM.options.han
%             % apply hanning window
%            P = P .* ( (1 - cos(2*pi*[1:Nf]'/(Nf + 1)))/2 )'; 
%         end
        
        try
            DCM.xY.csd{i}=P';
        catch
            DCM.xY.csd{i}=P;
        end   
end
    
%end

 
% place cross-spectral density in xY.y
%==========================================================================
DCM.xY.y    = spm_cond_units(DCM.xY.csd,'csd'); 
% 
% %KRISH HACK! If No normalisation is selected, do not normalise by power
   DoNormalise=1;
    try
        DoNormalise=DCM.options.DoNormalise;
        fprintf(1,'\nDoNormalise=%f\n',DoNormalise);
    catch
        DoNormalise=1;
    end

 if DoNormalise==0;
     DCM.xY.y    = krish_cond_unitsNONORMALISE(DCM.xY.csd,'csd'); 
     fprintf(1,'Not normalising to power....\n');
 end
 
 if DoNormalise==2;
     DCM.xY.y    = DCM.xY.csd; 
     fprintf(1,'Simple Spectral Copy....\n');
 end

 
 
DCM.xY.U    = DCM.M.U;
DCM.xY.code = condlabels(trial);
