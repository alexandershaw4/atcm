function timefreq = bert_singlechannel(InputData, cfg, freqvec, Basetimes)

%
% function timefreq = bert_singlechannel(InputData, [cfg], [freqvec], [Basetimes])
%
% Perform Hilbert time-frequency analysis of a single channel
%
% INPUTS 
%        InputData either
%        a) 1 D Matrix containing a single time-series to analyse 
%        b) 2 D matrix (Trials * Samples) of Data to Analyse
%        c) 3 D Matrix (Trials * Channels * Samples)
%        d) A fieldtrip Data structure (e.g. from makeSAMsensors)
%
%        cfg - a configuration structure with various options detailed below
%            - omit or use [] for defaults
%            - Remember if are using matrix style input data then cfg.fsample must be specified!  
%
%        freqvec - specify a vector e.g. 1:1:100 of start frequencies to analyse
%                - this is a shorthand way to specify cfg.freqvec   
%                - can be omitted or set as []
%
%        Basetimes - a 2 unit vector e.g [-1 0]
%                  - this is a shorthand way for entering cfg.base_start and cfg.base_end
%
% OUTPUTS
%
%        timefreq = structure containing various output fields 
%
% DATA SPECIFICATION OPTIONS
%
% cfg.fsample - the sample rate of the data - if the input data is just a matrix then this must be specified
%
%  
% cfg.sampletimes or (cfg.start_time and cfg.end_time)
%      - these options are for non fieldtrip data
%        cfg.sampletimes would be a vector containing the time values of each sample index therefore it would be as long as the number of samples in the data 
%        cfg.start_time and cfg.endtime e.g. cfg.start_time = -1 and cfg.endtime = 2 (a -1 s to 2 s trial) 
%      - If neither of these options are used then an arbitrary sampletimes vector will be created from time 0 taking into account the
%        the sampling frequency i.e . 0 :  1/ Fs : end 
%
%
% For nonfieldtripdata either Sampletimes or  
%     If fieldtrip data is ented this is used
%
% cfg.venumber = which sensor/virtual electrode to analyse for input type
%               a)  no selection is needed
%               b) no selection is needed
%               c) If a 3D matrix is entered then this selects the channel
%               d) If a fieldtrip structure is input this selects the appropriate channel
%               For c) and d) if novenumber is selected default = 1
%
% FREQUENCY OPTIONS
%
% cfg.bandwidth choose the bandwidth of your analysis - this defaults to 8 Hz
% 
% Specification of frequency can be done by specifying either
% 
% a) cfg.startf cfg.endf or cfg.fstep
% b) cfg.freqvec 
% c) cfg.freqs = the actual centre frequencies to be analysed (this is the equivalent cfg.freqvec + (bandwith / 2) 
% d) specify nothing and a cfg.freqvec of 1:1:100 will be used as default
%
%
% BASELINE OPTIONS
%
%  cfg.baseline = 'relchange' - the default but this should not be used for statistical purposes
%               = 'subtract'  - subtract the baseline epoch form the data 
%               = 'none'      - do nothing
%
%  cfg.base_start and cfg.base_end
%          Specify when to apply the baseline times e.g. cfg.base_start = -2 and cfg.base_end = 0
%          If not specified interval becomes [maxneg 0] or if no negative times just the first 100 samples
%
%
%
% STATISTICAL OPTIONS
%
%  cfg.ttest = 'yes' ; compute a t statistic spectrogram
%                       - the tgram/tgram_pval  are not corrected for multiple pcomparisons (use fdr for this perhaps)
%
% cfg.permutogram = 1000 ;construct a permutation p value gram with 1000 permutations. This returns
%               - pgram   - the pixelwise probability spectrogram           
%               - maxhist - the histogram of maximum permuted values
%               - minhist - the histrogram of minimum values
%               -cfg.permutogram and cfg.keeptrials are mutually exclusive
%     If you dont understand what these things mean you need to read Nichols and Holmes's paper       
%     To compute this the entire trials x time x frequency structure must be held in memory so you need to tihnk about memory issues! 
%     My playing around suggested that relative change baselines are very volatile so a subtract baseline is used for computation the histograms threfore
%     reflect the output (max/min values of a relative change histogram. However, by nature the permutations have more power/variance in the low frequencies and
%     the max/min hists probably simply reflect this variance. I would suggest that max/min distributions/omnibus test statstic are not appropriate to use
%     My suggestion is to compute your agram using relchange baseline so that it looks good      
%     The permutation testing will be done (using a subtract baseline) and then you can use
%      [p_masked, p_fdr] = fdr(1 - abs(timefreq.pgram), 0.05) % which will give an FDR corrected voxelwise significance mask, that you can apply to the relchange gram :-) 
%      finalgram = p_fdr .* timefreq.agram ;               
%               I havent directly implemented this because that way you can bodge your alpha level accordingly; well this neuroimaging after all :-)
%
%
% OSCDIAG OPTIONS
%
%    Oscdiag is a tool to test the reliabiltiy of (peak) frequency
%    estimation in a time-frequency spectrum (e.g.agram), e.g. returns
%    confidence intervals around the peak frequenncy. To run oscdiag on
%    your data set up an oscdiag substructure in cfg , that is cfg.oscdiag
%
%    You must specify
%    cfg.oscdiag.stimstart - Start time of time-frequency ROI e.g. cfg.oscdiag.stimstart= 0;
%    cfg.oscdiag.stimend   - End time of time-frequency   ROI e.g. cfg.oscdiag.stimend= 0.3;
%    cfg.oscdiag.stimfreqlow  - Lower frequency range to exmaine e.g. cfg.oscdiag.stimfreqlow = 30;
%    cfg.oscdiag.stimfreqhigh - Upper frequency range of interest e.g. cfg.oscdiag.stimfreqhigh = 100; 
%
%    You can specify as many time-frequency ROIs by making each of the 4
%    above variables vectors. Obviously the 4 vectors should be of even
%    length!
%
%    cfg.oscdiag.NPerms - The number of bootstrap iterations to run. This is optional and defaults to 5000    
%
%    cfg.oscdiag.CI - The confidence Interval boundaries you want e.g. cfg.oscdiag.CI = [0.025 0.975] (which is default) - i.e. 95 %
%    confidence interval
%
%    If oscdiag functionality is running too slow you may wish to consider pre-decimation of the data or better yet
%    use AllBerts_cluster to run it on the cluster (and select your channel of interest
%
% ADVANCED OPTIONS
% 
% cfg.filterorder - this defaults to 3 if you dont what this means then dont change it :-)
%
% cfg.predecimate e.g. cfg.predecimate = 2 ; Decimate the data by a factor of 2 prior to computation
%                - this will speed up computation time
%
% cfg.postdecimate e.g. cfg.postdecimate = 2 ; decimate the output by a factor of 2
%                - Computation time wont really be affected but this will save memory for big statistical gram operations!
%
% cfg.compute - 
%         'abs' -absolute value of the Hilbert trace (default) 
%         'pow' - Compute the power of the Hilbert trace
%         'phase' - COmpute the phase of the Hilbert trace     
%         'complex' keep the complex elements
%         'plv' - compute the PLV across
%
% cfg.central
%        'mean' - default (of the agram)
%        'median'   %median only makes sense for abs/pow/phase
%
% cfg.kill_evagram
%        1 = dont save evagram
%        0 = save it (if not specfied will be saved
%
% cfg.keeptrials - keep the individual trials x time x frequency structure
%                - you need to think about how much memory this is going to consume
%                - consider using pre or post decimation
%                - data is returned unbaselined
%                -cfg.premutogram and cfg.keeptrials are mutually exclusive 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Some example usages
% 1. timefreq = bert_singlechannel(FieldTripData,[],10:0.5:100,[-2 0])
%      This does a timefrequency data of FieldTripData (from makeSAMsensors).
%      The agram is computed from 10 to 100 Hz in 0.5 Hz steps.
%      A default relative change baseline will be used using -2 to 0 as the baseline times
%
% 2. cfg.fsample = 600;
%    cfg.ttest = 'yes'
%    cfg.compute = 'pow'
%    cfg.central = 'median'
%    cfg.start_time = -2
%    cfg.endtime = 2;  
%    cfg.venumber = 3;
%    timefreq = bert_singlechannel(MatrixData,cfg,10:0.5:100,[-2 0])
%       This does time frequency analysis of the Data held in  MatrixData
%       Since it is Matrix Data the sample rate must be specfied
%       We also specify the start/endtimes so it can generate sampletimes for us
%       The 3rd electrode in the structure is analysed (venumber)
%       THe agram will contain a median spectrogram of the power and there will also be a ttest spectrogram for us
%
% 3. cfg.compute = 'plv';
%    cfg.baseline = 'none'
%    timefreq = bert_singlechannel(FieldTripData,cfg,10:0.5:100)
%       Similar to Example 1 but here the PLV (a measure of phase consistency across trials) will be calculated
%       At least > 20 trials are needed for the PLV computation to even slightly work....
%
% 4. cfg.permutogram = 1000;
%    cfg.postdecimate = 4;
%    cfg.baseline = 'subtract'
%    cfg.ttest = 'yes'
%    timefreq = bert_singlechannel(FieldTripData,cfg,10:0.5:100,[-2 0])
%       Compute the time-fequency spectrogram using a subtracted baseline
%       Compute a ztest spectrogram
%       Also compute a permutation spectrogram
%       All the outputs will be downsampled by a factor of 4 after computation
%
% Some post-processing tips
%   1. Use GetTimeFreqData to extract numbers from the time-frequency data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Testing
%   load CPA1_Pre_V1
%   InputData = CPA1_Pre_V1_Good
% %  cfg = [];
% % 
% 
% %load InputData
% nargin  = 2;
% cfg = [];
% cfg.fsample = 1200;
% cfg.bandwidth = 8;    
% cfg.startf = 1;
% cfg.endf = 80;
% cfg.fstep = 2;
% cfg.compute = 'abs';
% cfg.central = 'mean';
% cfg.baseline = 'subtract';
% cfg.basestart = -2
% cfg.postdecimate = 4;
% cfg.permutogram = 500;
%cfg.ttest = 'yes'
%cfg.freqvec = 1 : 1: 50
%cfg.freqs = 5:1:104

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if checkcheck ~= 1234567890, error('Bad cheque'); end

if nargin == 0
    help bert_singlechannel
    return
end


timefreq = [];

if isstruct(InputData) == 1
    IsFieldTrip = 1;
else
    IsFieldTrip = 0;
end



if nargin == 1
   cfg = []; 
end

%Sort out the sample rate
if IsFieldTrip == 1
    Fs = InputData.fsample;
else
   if isfield(cfg, 'fsample') == 1
       Fs = cfg.fsample;
   else
       error('For non fieldtrip data a sample rate must be specified (e.g. cfg.sample = 1200)');
   end
end

%Sort out the Channel specified to analyse
if isfield(cfg, 'venumber') == 1
    venumber = cfg.venumber;
    fprintf('Analysing channel %d in the data\n', venumber);
else
    venumber = 1;
    %fprintf('No venumber was specfied so analysing the first channel in the data\n');
end

%Sort the Data out
% Create "vedata" - a 2d matrix that contains the data to analyse - epochs x samples
if IsFieldTrip == 1 %Extract the data out of the fieldtrip structure and put into vedata ready for analysis
   NEpochs = length(InputData.trial);
   NSamples = size(InputData.trial{1}, 2);
   vedata = zeros([NEpochs, NSamples]);
   for m = 1 : NEpochs
      vedata(m,:) = InputData.trial{m}(venumber,:);
   end
else
    if ndims(InputData) == 1
        vedata = InputData;
        NEpochs = 1;
        NSamples = length(vedata);
    elseif ndims(InputData) == 2
        vedata = InputData;
        [NEpochs NSamples] = size(vedata);      
    elseif ndims(InputData) == 3
        vedata = squeeze(InputData(:,venumber,:));
        [NEpochs NSamples] = size(vedata);
    end
end

%Sort out the outputtime base 'sampletimes'
if IsFieldTrip == 1 
    if iscell(InputData.time)
       sampletimes=cell2mat(InputData.time(1));
    else
       sampletimes = InputData.time;   
    end
else 
    if isfield(cfg, 'sampletimes')
        sampletimes = cfg.sampletimes;
    elseif isfield(cfg, 'start_time') && isfield(cfg, 'end_time')
        sampletimes = linspace(cfg.start_time, cfg.end_time, NSamples);
    else
        fprintf('An arbitrary output time vector is being generated\n');
        sampletimes = 0: 1/ Fs : ((NSamples / Fs)- (1/Fs));
    end
end

%Sort Out the Baselining
if isfield(cfg, 'baseline') == 0
    baseline = 'relchange';
    %fprintf('Using a relative baseline as default\n');
else
    baseline = cfg.baseline; 
end

if (nargin > 3)
    cfg.base_start = Basetimes(1);
    cfg.base_end = Basetimes(2);
end
    
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Frequency paramters

if isfield(cfg, 'filterorder') == 0
   FilterOrder = 3;
else
   FilterOrder = cfg.filterorder;
end

if isfield(cfg, 'bandwidth') == 0
   bandwidth = 8;
else
   bandwidth = cfg.bandwidth;
end

if (nargin > 2) && (isempty(freqvec) == 1) || (nargin <= 2) %If freqvec is given and not specified or omitted
    if isfield(cfg, 'startf') == 1 && isfield(cfg, 'endf') == 1 && isfield(cfg, 'fstep') == 1
       freqvec = cfg.startf : cfg.fstep : cfg.endf ;
    elseif isfield(cfg, 'freqvec') == 1
       freqvec = cfg.freqvec;
    elseif isfield(cfg, 'freqs') == 1
       freqvec = cfg.freqs - (bandwidth / 2);
    else
       freqvec = 1 : 1 : 100;
    end
end

%Now create the frequency information for the analysis
NFreqs = length(freqvec);
for m = 1 : NFreqs
   hpHz(m)= freqvec(m);
   lpHz(m)= hpHz(m) + bandwidth;
   freqs(m) =(hpHz(m)+lpHz(m))/2; %freqs contains the centre
end

%What statistic are we goint to compute
if isfield(cfg, 'compute') == 0
    compute = 'abs';
    %fprintf('Computing the absolute value of the Hilbert trace by default\n');
else
    compute = cfg.compute;
end
if isfield(cfg, 'central') == 0
    central = 'mean';
    %fprintf('Computing the mean value of computation variable by default\n');
else
    central = cfg.central;
end
if isfield(cfg, 'ttest') == 1
    DoTTest = 1;
    fprintf('Will compute a t statistic spectrogram\n');
else
    DoTTest = 0;
end

if isfield(cfg, 'keeptrials') && isfield(cfg, 'permutogram')
   error('You cannot specify both keeptrials and permutogram - sorry'); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Start the actual computations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%If requested here we predecimate the data and adjust variables accordingly ....
if isfield(cfg, 'predecimate')
   fprintf('Predecimating the data by a factor of %d\n', cfg.predecimate);
   for m = 1 : NEpochs,
      vedatanew(m,:)=decimate(vedata(m,:),cfg.predecimate); %Decimate to avoid alias
   end
   vedata = vedatanew;
   sampletimes=downsample(sampletimes,cfg.predecimate); %Just downsample this
   Fs=Fs/cfg.predecimate;
   [NEpochs NSamples] = size(vedata);
end

%Now the predecimation is done we can work out when the baselining will be in terms of samples
if strcmp(baseline, 'none') == 1
else
    if isfield(cfg, 'base_start') && isfield(cfg, 'base_end') ==1
        base_start_time =  cfg.base_start;
        base_end_time =  cfg.base_end;

        base_start = findthenearest(base_start_time,sampletimes,1);
        base_end = findthenearest(base_end_time,sampletimes,0);
    else
        if sampletimes(1) < 0 && sampletimes(end) > 0 %If it is a stand time with a negative baseline e.g. -2 to 0
           base_start_time = sampletimes(1);
           base_start = 1;
           base_end_time = 0;
           base_end = findthenearest(0,sampletimes,-1);
        else
           base_start_time = sampletimes(1);
           base_start = 1;
           base_end_time = sampletimes(100);
           base_end = 100;
        end
       %fprintf('Using the baseline period of samples %d to %d\n',base_start, base_end );
    end
end
%%%%%%%%%%%%%%%%

if isfield(cfg, 'postdecimate')
   PostDecimate = cfg.postdecimate;
   SamplesToKeep = 1:PostDecimate:length(sampletimes); %These are used for keeping the 3D trial info
   NSamplesToKeep = length(SamplesToKeep);             %These are used for keeping the 3D trial info
else
   PostDecimate = 1; 
   SamplesToKeep = 1 : 1 : length(sampletimes);
   NSamplesToKeep = NSamples;
end

%First we should pad the data to avoid edge effects (the data we will analyse is called "extended_data"
extension = 0.2;
buflength = round(extension*NSamples);
NExtended_Samples = NSamples + (2 * buflength);
extended_data = zeros([NEpochs NExtended_Samples]);
extended_data(:, buflength+1: buflength+ NSamples)=vedata;  %Move the vedata into the middle
extended_data(:,1:buflength) = fliplr(vedata(:,2:buflength+1)); %Buffer the start by mirroring the start data
extended_data(:,buflength+ NSamples+1:end) = fliplr(vedata(:,NSamples - buflength :end-1)); %Buffer the end by mirroring the end data
                                     %After computation we only want to keep buflength+ 1 : buflength + NSamples

%%%%%%%%%%%%%%%Now lets calculate the 'Evoked'
if NEpochs > 1
   evoked = mean(extended_data);
end

%%%
if isfield(cfg, 'permutogram') || isfield(cfg, 'keeptrials') || isfield(cfg, 'oscdiag')
   KeepTrials = 1;
   fprintf('Attempting to get memory for Large 3 D matrix - hold onto your hats......');
   timefreq.Allgrams = zeros([NEpochs, NFreqs, NSamplesToKeep]); 
   fprintf('Done\n');
else
    KeepTrials = 0;
end
    

temp_gram = zeros([NEpochs NSamples]);

evagram = zeros([NFreqs NSamples]);
Outgram = zeros([NFreqs NSamples]);

if  DoTTest == 1
   tgram = zeros([NFreqs NSamples]);
   tgram_pval = zeros([NFreqs NSamples]);
   
end
    

%The main loop is frequency
for m = 1 : NFreqs
    %if m > 1; fprintf(repmat('\b',size(str))); end
    %str = sprintf('Computing frequency %d of %d \n', m, NFreqs);
    %fprintf(str);
    
    %fprintf('Computing frequency %d of %d \n', m, NFreqs);
    %fprintf('[computing %s statistic]\n',compute);
    %Design the filter coefficients 
    [B, A] = atcm.fun.mk_filter(Fs,hpHz(m), lpHz(m), FilterOrder);
       for n = 1 : NEpochs
             filtve=filtfilt(B,A,extended_data(n,:));
             trace= hilbert(filtve);
             temp_gram(n,:) = trace(buflength+ 1 : buflength + NSamples);  %tempgram holds all the complex numbers for this frequency (Epochs x Samples
       end
       %Now compute the relevant measure from the complex numbers
       if strcmp(compute, 'abs') ==1
          temp_gram = abs(temp_gram);
       elseif strcmp(compute, 'logabs') ==1 %This is Suresh fudging about
          temp_gram = log10(abs(temp_gram));
       elseif strcmp(compute, 'pow') ==1   
          temp_gram = abs(temp_gram).^2;
       elseif strcmp(compute, 'phase') ==1 || strcmp(compute, 'plv') ==1   
          temp_gram = angle(temp_gram);
       end
       %If we need to keep the numbers store them
       if KeepTrials == 1
           timefreq.Allgrams(:,m,:) = temp_gram(:,SamplesToKeep); %This implements a quick post-decimation of the data (dont really want to hang onto all of it)
       end
       %Now compute our measures of central tendency
       if NEpochs > 1
           if strcmp(compute, 'plv') ==1 %For PLV calculation there can be no mean/median 
              Outgram(m,:) = (abs(sum(exp(temp_gram * 1i)))) ./ NEpochs;
           else
              if (strcmp(compute, 'abs') || strcmp(compute, 'pow') || strcmp(compute, 'complex')|| strcmp(compute, 'logabs') ) && strcmp(central, 'mean') %Cant say I am 100% sure what the mean complex gives you!
                 Outgram(m,:) = mean(temp_gram);
              elseif (strcmp(compute, 'abs') || strcmp(compute, 'pow')|| strcmp(compute, 'complex')) && strcmp(central, 'median')   
                 Outgram(m,:) = median(temp_gram);
              elseif (strcmp(compute, 'phase') && strcmp(central, 'mean'))
                 Outgram(m,:) = circ_mean(temp_gram);
              elseif (strcmp(compute, 'phase') && strcmp(central, 'median'))
                 Outgram(m,:) = circ_median(temp_gram);
              end
           end
       
           %Compute the evagram for this frequency
           filtve=filtfilt(B,A,evoked);
           trace= hilbert(filtve);
           if strcmp(compute, 'pow') == 1
              evagram(m,:) = abs(trace(buflength+ 1 : buflength + NSamples)) .^ 2;
           elseif strcmp(compute, 'phase') == 1
              evagram(m,:) = angle(trace(buflength+ 1 : buflength + NSamples));
           else %If not power or phase just compute the abs for the evagram
               evagram(m,:) = abs(trace(buflength+ 1 : buflength + NSamples));
           end
       else   %FOr single trial data
           Outgram(m,:) = temp_gram;
       end

       %We can compute t statistic grams "tgrams" by frequency. However permutation testing can only be done when we have the entire 3 DImensional time-fequency array....
       if DoTTest == 1 %This works for power and abs only!
          % - First subtract baseline from 
          tempgram_base = Baseline2D(temp_gram, base_start, base_end); %THe Baseline needs to be applied else its gonna look spazzy :-)
          [H,P,CI,STATS] = ttest(tempgram_base);
          tgram(m,:) = STATS.tstat;
          tgram_pval(m,:) = 1 -P; % 1- p value images plot better :-)
          %Create t test image
          
       end
       
end

%Now we are finished with evoked shrink it
if NEpochs > 1
   evoked = evoked(buflength+ 1 : buflength + NSamples);
end

%Now that we have computed the Spectrogram/Evagram Baseline It
if strcmp(baseline, 'relchange') ==1 
   Outgram = atcm.fun.Baseline2D(Outgram, base_start,base_end, 0);
   if NEpochs > 1
      evagram = atcm.fun.Baseline2D(evagram, base_start,base_end, 0);
   end
elseif strcmp(baseline, 'subtract') ==1 
   Outgram = atcm.fun.Baseline2D(Outgram, base_start,base_end);
   if NEpochs > 1
      evagram = atcm.fun.Baseline2D(evagram, base_start,base_end);
   end
end

%Shrink all the variables down just to what we are going to keep decided by the PostDecimate factor
if PostDecimate > 1
    sampletimes=downsample(sampletimes,PostDecimate);
    evoked = decimate(evoked, PostDecimate);
    Fs = Fs / PostDecimate;
    for m = 1 : NFreqs,
       Outgram_temp(m,:) = decimate(Outgram(m,:), PostDecimate);
       evagram_temp(m,:) = decimate(evagram(m,:), PostDecimate);
    end
    Outgram = Outgram_temp;
    evagram = evagram_temp;
    if DoTTest == 1
        for m = 1 : NFreqs,
            tgram_temp(m,:)   = decimate(tgram(m,:), PostDecimate);
            tgram_pval_temp(m,:)   = decimate(tgram_pval(m,:), PostDecimate);
        end    
        tgram = tgram_temp;
        tgram_pval = tgram_pval_temp;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Do Oscdiag bootstrapping here
if isfield(cfg, 'oscdiag')
    %Get the Number to Permute
    if isfield(cfg.oscdiag, 'NPerms')
        osc_reps = cfg.oscdiag.NPerms;
    else
        osc_reps = 5000;
    end
    %Get the number of time-frequency bins top analyse and error check
    osc_tfwindows = length(cfg.oscdiag.stimstart);
    fprintf('There are %d oscdiag time-frequency windows to compute\n', osc_tfwindows);
    if (length(cfg.oscdiag.stimend) ~= osc_tfwindows) || (length(cfg.oscdiag.stimfreqlow) ~= osc_tfwindows) || (length(cfg.oscdiag.stimfreqhigh) ~= osc_tfwindows)
        error('oscdiag.stimstart oscdiag.stimend oscdiag.stimfreqlow oscdiag.stimfreqhigh must all be the same length');
    end
    %Get the indexes of the time-frequency bins for array referencing
    for m = 1 : osc_tfwindows      
        osc_stimstart(m) =findthenearest(cfg.oscdiag.stimstart(m),sampletimes,1);
        osc_stimend(m) =findthenearest(cfg.oscdiag.stimend(m),sampletimes,-1);
        osc_stimfreqlow(m) =findthenearest(cfg.oscdiag.stimfreqlow(m),freqs,1);
        osc_stimfreqhigh(m) =findthenearest(cfg.oscdiag.stimfreqhigh(m),freqs,-1);
        oscdiag{m}.freqs = freqs( osc_stimfreqlow(m): osc_stimfreqhigh(m));
    end
    %We can extract the obtained mean value from Outgram (this has already had baseline correction performed)
    for m = 1 : osc_tfwindows
         oscdiag{m}.meanspectrum = mean(Outgram(osc_stimfreqlow(m):osc_stimfreqhigh(m),osc_stimstart(m):osc_stimend(m)),2);
         oscdiag{m}.newfreqs     = oscdiag{m}.freqs(1):0.01: oscdiag{m}.freqs(end);
         oscdiag{m}.meanspectrum = spline(oscdiag{m}.freqs,oscdiag{m}.meanspectrum, oscdiag{m}.newfreqs);
         
         [oscdiag{m}.MaxValue,MaxIndex] = max(oscdiag{m}.meanspectrum);
         [oscdiag{m}.MinValue,MinIndex] = min(oscdiag{m}.meanspectrum);
         oscdiag{m}.MaxFreq = oscdiag{m}.newfreqs(MaxIndex);
         oscdiag{m}.MinFreq = oscdiag{m}.newfreqs(MinIndex);
    end     
    %First allocate the memory for the outputvariables
    for m = 1 : osc_tfwindows
        oscdiag{m}.ArrayMaxFreq=(1:osc_reps);
        oscdiag{m}.ArrayMaxValue=(1:osc_reps);
        oscdiag{m}.ArrayMinFreq=(1:osc_reps);
        oscdiag{m}.ArrayMinValue=(1:osc_reps);     
    end
    
    osc_indices = squeeze(1+fix(rand(osc_reps,NEpochs).*(NEpochs))); % generates random sample indices
    %Now we need do the bootstrapping
    for n = 1 : osc_reps
        %Compute the permuted mean and apply baseline correction
        fprintf('oscdiag - running permutation %d of %d\n', n, osc_reps);
        osc_permgram = squeeze(mean(timefreq.Allgrams(osc_indices(n,:),:,:),1));
        if strcmp(baseline, 'relchange') ==1 
             osc_permgram = Baseline2D(osc_permgram, base_start,base_end, 0);
        elseif strcmp(baseline, 'subtract') ==1 
             osc_permgram = Baseline2D(osc_permgram, base_start,base_end);
        end
        %Grab the Variables we need 
        for m = 1 : osc_tfwindows
             osc_meanspectrum = mean(osc_permgram(osc_stimfreqlow(m):osc_stimfreqhigh(m),osc_stimstart(m):osc_stimend(m)),2);
             osc_meanspectrum = spline(oscdiag{m}.freqs,osc_meanspectrum, oscdiag{m}.newfreqs);
             [oscdiag{m}.ArrayMaxValue(n),MaxIndex] = max(osc_meanspectrum);
             [oscdiag{m}.ArrayMinValue(n),MinIndex] = min(osc_meanspectrum);
             oscdiag{m}.ArrayMaxFreq(n) = oscdiag{m}.newfreqs(MaxIndex);
             oscdiag{m}.ArrayMinFreq(n) = oscdiag{m}.newfreqs(MinIndex);
        end    
    end
    
    %Now calculate confidence intervals and clean-up
    if isfield(cfg.oscdiag, 'CI')
        osc_lowCI = cfg.oscdiag.CI(1);
        osc_highCI = cfg.oscdiag.CI(2);
    else
        osc_lowCI = 0.025;
        osc_highCI = 0.975;
        fprintf('Default confidence intervals of 0.025 to 0.975 will be generated\n');
    end
    osc_lowthresh = round(osc_lowCI * osc_reps);
    osc_highthresh = round(osc_highCI * osc_reps);
    for m = 1 : osc_tfwindows
       osc_temp = sort(oscdiag{m}.ArrayMaxValue);
       oscdiag{m}.MaxValueCI = [osc_temp(osc_lowthresh) osc_temp(osc_highthresh)];
       osc_temp = sort(oscdiag{m}.ArrayMinValue);
       oscdiag{m}.MinValueCI = [osc_temp(osc_lowthresh) osc_temp(osc_highthresh)];
       osc_temp = sort(oscdiag{m}.ArrayMaxFreq);
       oscdiag{m}.MaxFreqCI = [osc_temp(osc_lowthresh) osc_temp(osc_highthresh)];
       osc_temp = sort(oscdiag{m}.ArrayMinFreq);
       oscdiag{m}.MinFreqCI = [osc_temp(osc_lowthresh) osc_temp(osc_highthresh)];
       oscdiag{m}.sampletimes = sampletimes; 
    end     
    timefreq.oscdiag = oscdiag;
    clear osc*
    fprintf('Finished running oscdiag options\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Permutation Testing Code - this only works for abs and power

%IMPORTANT note - the permutation testing code shuffles the big 3D agram
%and doesnt keep track of where it is shuffled to - so it needs to be
%discarded!!
%
if isfield(cfg, 'permutogram')
    fprintf('Begining permutation testing\n');
    NPerms = cfg.permutogram;
    if strcmp(compute, 'abs') == 1
        timefreq.Allgrams = abs(timefreq.Allgrams);
    elseif strcmp(compute, 'pow') == 1    
        timefreq.Allgrams = abs(timefreq.Allgrams) .^ 2;
    end
    %The Baseline samples need to be adjusted because now we are post decimation
    base_start = findthenearest(base_start_time,sampletimes,1);
    base_end = findthenearest(base_end_time,sampletimes,-1);
    
    %We need a thispermgram variable to hold onto the results of this permutation
    timefreq.maxhist = zeros([1 NPerms]);
    timefreq.minhist = zeros([1 NPerms]);
    timefreq.pgram = zeros([NFreqs NSamplesToKeep]);
    countgram = zeros([NFreqs NSamplesToKeep]);
    
    %Lets recompute the obtained spectrogram in case a relchange baseline was used
    if strcmp(central, 'mean') == 1
       obtainedgram= squeeze(mean(timefreq.Allgrams));
    elseif strcmp(central, 'median') == 1
       obtainedgram= squeeze(median(timefreq.Allgrams));
    end
    obtainedgram = Baseline2D(obtainedgram, base_start,base_end);
    
       
    for k = 1 : NPerms
        fprintf('Computing permutation %d of %d \n', k, NPerms);
        permvec = (round(rand([1 NEpochs]))*2)-1; %compute a random permutation vector of -1 and 1 s
        %Permute the sign of the trials
        for m = 1 : NEpochs
           timefreq.Allgrams(m,:,:) = timefreq.Allgrams(m,:,:) * permvec(m);  
        end
        %Compute the mean or median of the variable
        if strcmp(central, 'mean') == 1
           thispermgram = squeeze(mean(timefreq.Allgrams));  
        else
           thispermgram = squeeze(median(timefreq.Allgrams));  
        end
        %Then we need to baseline the matrix - goning to force subtract baselining
        thispermgram = Baseline2D(thispermgram, base_start,base_end);
        %Now get the max and mins for this permutation
        timefreq.maxhist(k) = max(max(thispermgram)); 
        timefreq.minhist(k) = min(min(thispermgram));
        %There has got be a faster way of doing the following !!        
        for m = 1 : NFreqs
            for n = 1 : NSamplesToKeep
                if (obtainedgram(m,n) > 0) && (thispermgram(m,n) > obtainedgram(m,n))    %If its a positive time-f value and bigger than obtained count it
                     countgram(m,n) = countgram(m,n) + 0.5; %The 0.5 makes it two-tailed
                elseif (obtainedgram(m,n) < 0) && (thispermgram(m,n) < obtainedgram(m,n))  %If its a negative time-f value but more negative count it
                     countgram(m,n) = countgram(m,n) - 0.5;
                end
            end
        end
    end%ENd of the permutations
    %Create the pgram from the countgram
    for m = 1 : NFreqs
        for n = 1 : NSamplesToKeep
            if obtainedgram(m,n) > 0
                timefreq.pgram(m,n) = (NPerms - countgram(m,n)) /NPerms;
            else
                timefreq.pgram(m,n) = (-NPerms - countgram(m,n)) /NPerms;
            end
        end
    end
     
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(cfg, 'keeptrials') && isfield(timefreq, 'Allgrams')
       timefreq.Allgrams = [] ;  %Kill the 3D gram now we are done with it
end

%%%%%%Gather up the ouputs
if DoTTest ==1
   timefreq.tgram = single(tgram);
   timefreq.tgram_pval = single(tgram_pval);
end
        
timefreq.agram = single(Outgram);
timefreq.sampletimes = single(sampletimes);
timefreq.freqs = single(freqs);
timefreq.Fs = Fs;
timefreq.nepochs = NEpochs;
timefreq.cfg = cfg; %Lets store this for user reference!!

if NEpochs > 1
    if ~isfield(cfg, 'kill_evagram') || cfg.kill_evagram ~= 1
        timefreq.evagram = single(evagram);      
    end
     timefreq.evoked = single(evoked);
end


%fprintf('bert_singlechannel complete\n');
%All done!

%%%%%%%%%%%%%%%%%%%%%%


