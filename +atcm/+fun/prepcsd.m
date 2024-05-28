function DCM = prepcsd(DCM)
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
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_ssr_data.m 4096 2010-10-22 19:40:34Z karl $
 
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

% alex - allow [] to specify all
if isempty(DCM.options.trials)
    trial = 1;
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
    if ~isempty(DCM.options.trials)
        c = D.indtrial(condlabels(trial(i)), 'GOOD');
    else
        c = 1:size(D,3); 
    end
    
    if isfield(DCM.options,'TrialIndices')
        fprintf('Using only specified trial range...\n');
        c = c(DCM.options.TrialIndices);
    end
    
    
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
        
        % log the condition and trial series
        %series{i}(j,:,:) = Y;
        
        if DOBASE==1;
                Ybase=full(double(D(Ic,baseIt,c(j))'*DCM.M.U));
        end
            
        try
            DoButterband=DCM.options.UseButterband;
            %Y=butterband(Y,DoButterband(1),DoButterband(2),Fs);
            
            for ch = 1:size(Y,2)
                %[dY,in] = atcm.fun.padtimeseries(Y(:,ch));
                %dY      = atcm.fun.bandpassfilter(dY',1/DCM.xY.dt,DoButterband)';
                %xY      = dY(in);
                %Y(:,ch) = xY;
                Y(:,ch) = atcm.fun.bandpassfilter(Y(:,ch),1/DCM.xY.dt,DoButterband)';
            end 
            
            if DOBASE
                for ch = 1:size(Ybase,2)
                    %[dYb,in] = atcm.fun.padtimeseries(Ybase(:,ch));
                    %dYb      = atcm.fun.bandpassfilter(dYb',1/DCM.xY.dt,DoButterband)';
                    %xY       = dYb(in);
                    %Ybase(:,ch)  = xY;
                    Ybase(:,ch) = atcm.fun.bandpassfilter(Ybase(:,ch),1/DCM.xY.dt,DoButterband)';
                end
            end
            
        end
        
%         % baseline
%         if DOBASE
%             %fprintf('Baselining\n');
%             BaseData = repmat( mean(Ybase,1), [1,ch] );
%             Ymod = (Y - BaseData);%./BaseData;
%         else
%             Ymod = Y;
%         end
        Ymod = Y;
        series{i}(j,:,:) = Y;
        
        UseWelch=0;
        try 
            UseWelch=DCM.options.UseWelch;
        catch
            UseWelch=0;
        end
        if(UseWelch==1010)
           
            if isfield(DCM.options,'vmd') && DCM.options.vmd
                [Ymod] = vmd(real(Ymod),'NumIMFs',5);
                %fprintf('Using VMD\n');
                
                [uu,~,~] = spm_svd(cov(Ymod'));                                        
                pcc = uu'*Ymod;
                nn = min(12,size(pcc,1));%thr;
                pcc = pcc(1:nn,:);
                
                % autoregressive spectral method for VMD components
                for ii = 1:nn; Pfc(ii,:) = pyulear(pcc(ii,:),2,DCM.xY.Hz,1/DCM.xY.dt); end
                
                Pf(:,1,1) = sum(Pfc,1);
                
                 Pfull(j,:,:,:) = full(Pf);
                 
            else
            
            
                % Simple fft(x)(alex)
                %---------------------------------------------------
                FFTSmooth = 0;
                try FFTSmooth = DCM.options.FFTSmooth; end
                if UseWelch == 1010;

                    % Alex fft
                    SpecFun = DCM.options.SpecFun;
                    %[Pf, F] = atcm.fun.AfftSmooth(Ymod', 1/DCM.xY.dt, DCM.xY.Hz) ;

                    if strcmp(char(SpecFun),'atcm.fun.AfftSmooth') && isfield(DCM.options,'FFTbins');
                        [Pf, F] = SpecFun(Ymod', 1/DCM.xY.dt, DCM.xY.Hz,DCM.options.FFTbins) ;
                    else
                        [Pf, F] = SpecFun(Ymod',  1/DCM.xY.dt, DCM.xY.Hz) ;

                        % use the same regression based routine as the TCM
                        % model
                        if isfield(DCM.options,'RegressionDCM') && DCM.options.RegressionFFT
                            % if j == 1
                            %     [~,DFT,f] = adft(Ymod,1/DCM.xY.dt,DCM.xY.Hz);
                            % end
                            % 
                            % b = DFT\Ymod;
                            % 
                            % %[u,s,v] = svd(b);
                            % 
                            % %b = atcm.fun.Afft(ys,1/dt,w);
                            % 
                            % PpfR = atcm.fun.agauss_smooth(abs(real(b)),1);%u(:,1)'*b);
                            % PpfI = atcm.fun.agauss_smooth(abs(imag(b)),1);
                            % 
                            % Pf = PpfR + sqrt(-1)*PpfI;



                            if j == 1
                                DFT = lsqspectrum(Ymod,1/DCM.xY.dt,DCM.xY.Hz);
                            end

                            Pf = abs( (DFT)'\(ys)' );




                        end


                        if j == 1
                            FMAT   = atcm.fun.asinespectrum(DCM.xY.Hz,D.time(It)*1000);
                        end

                        %b   = atcm.fun.lsqnonneg(F',Ymod);
                        %b = abs(F'\Ymod);
                        %Pf  = atcm.fun.agauss_smooth(b,1);

                        %[Pf,F] = atcm.fun.afftg(Ymod,DCM.xY.dt,DCM.xY.Hz);

                        %[Pf,F] = atcm.fun.fftgc(Ymod',DCM.xY.dt,DCM.xY.Hz,10);

                        %F = DCM.xY.Hz;
                        %Pf = atcm.fun.tfdecomp(Ymod,DCM.xY.dt,F,2,3,@max);
                        %Pf = atcm.fun.approxlinfitgaussian(Pf,[],[],4);
                        
                        % if isfield(DCM.options,'fooof') && DCM.options.fooof
                        %     if j == Nt;fprintf('Fitting FOOOF %d times\n',j);end
                        %     Pf = atcm.fun.component_spectrum(F,Pf,12);
                        % end
                        
                        if size(Pf,1) == 1;
                            Pf = Pf';
                        end
                    end
%                     if FFTSmooth > 0;
%                         for nchanx = 1:size(Pf,2)
%                             for nchany = 1:size(Pf,3)
%                                 Pfxy = Pf(:,nchanx,nchany);
%                                 [dPf,in] = atcm.fun.padtimeseries(Pfxy);
%                                 dPf = atcm.fun.HighResMeanFilt(dPf,1,FFTSmooth);
%                                 Pfxy  = dPf(in);
%                                 Pf(:,nchanx,nchany) = Pfxy;
%                             end
%                         end
%                     end
                    % if isfield(DCM.options,'envelope') && ~isempty(DCM.options.envelope) && DCM.options.envelope>0
                    %    for nchanx = 1:size(Pf,2)
                    %        for nchany = 1:size(Pf,3)
                    %            Pf(:,nchanx,nchany) = atcm.fun.aenv(Pf(:,nchanx,nchany),DCM.options.envelope,[],1);
                    %            %Pf(:,nchanx,nchany) = atcm.fun.aenvelope(Pf(:,nchanx,nchany),DCM.options.envelope);
                    %        end
                    %    end
                    % end

                    Pfull(j,:,:,:) = full(Pf);
                end
            
            end % end vmd
            
            % the base period spectrum
            if DOBASE
                
                [Pfb, F] = atcm.fun.Afft(Ybase',  1/DCM.xY.dt, DCM.xY.Hz) ;
                if FFTSmooth > 0;
                    [dPfb,in] = atcm.fun.padtimeseries(Pfb);
                    dPfb = atcm.fun.HighResMeanFilt(dPfb,1,FFTSmooth);
                    Pfb  = dPfb(in);
                end
                Pfullbase(j,:,:,:) = full(Pfb);
                
            end
            
            
        end
    end
    
    
    % average trials for this condition
    if (UseWelch==1010)
        if isfield(DCM.options,'BeRobust') && DCM.options.BeRobust
            fprintf('Robust fitting\n');
            [mnewspectra,fq,unc,~,~,FitPar ]=atcm.fun.RobustSpectraFit(F,Pfull,2);
            mnewspectra=exp(mnewspectra);
            mnewspectra=mnewspectra-min(mnewspectra);
            Pfull = mnewspectra;

            %Pfull = spm_robust_average(Pfull);

            DCM.xY.trial_spectra = unc;
            
            if DOBASE
                Pfullbase = squeeze(spm_robust_average(Pfullbase));
            end
            
        else
            DCM.xY.trial_spectra = Pfull;
            %Pfull = squeeze(spm_robust_average(Pfull));
            
            if size(Pfull,1) > 1 && Nm==1
                % retain first eidenmode
                m        = 1;
                [u s v]  = spm_svd(Pfull',1);
                Pfull    = u(:,m)*s(m,m)*mean(v(:,m));
                Pfull    = full(Pfull)';
            elseif size(Pfull,1) > 1 && Nm>1
                Pfull = squeeze(spm_robust_average(Pfull));
            else
                Pfull = squeeze(Pfull);
            end
            
            if DOBASE
                Pfullbase = squeeze(spm_robust_average(Pfullbase));
            end
            
        end
        
        
        
        if FFTSmooth > 0;
            if Nm == 1
                Pfull = Pfull';
            end
            for nchanx = 1:size(Pf,2)
                for nchany = 1:size(Pf,3)
                    Pfxy = Pfull(:,nchanx,nchany);
                    [dPf,in] = atcm.fun.padtimeseries(Pfxy);
                    dPf = atcm.fun.HighResMeanFilt(dPf,1,FFTSmooth);
                    Pfxy  = dPf(in);
                    Pfull(:,nchanx,nchany) = Pfxy;
                end
            end
            if Nm == 1
                Pfull = Pfull';
            end
        end
        
        
 

        if DOBASE==1,
            %Pfullbase=squeeze(spm_robust_average(Pfullbase));
            Pfull = Pfull - Pfullbase;
        end
        
        P = Pfull;
        
        P(find(isinf(P)))=0;
        P(find(isnan(P)))=0;
        %P(find(P<=0))=0;
        P=P-min(P);
        
        if isfield(DCM.options,'han') && DCM.options.han
            % apply hanning window
            if Nc==1
                P = P .* rescale(kaiser(Nf,2.5),.01,1)'.^.2;
            else
                wind = rescale(kaiser(Nf,2.5),.01,1)';
                P = P .* repmat(wind(:),[1 Nc,Nc]).^.2;
            end
        end
        
        try
            DCM.xY.csd{i}=P';
        catch
            DCM.xY.csd{i}=P;
        end   
    end
end

 
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

% output the series
DCM.xY.series = series;
 
DCM.xY.U    = DCM.M.U;
DCM.xY.code = condlabels(trial);

try
    DCM.xY.FMAT = FMAT;
end
