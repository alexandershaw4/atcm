function [varargout] = vmd(x,varargin)
%VMD Variational mode decomposition
%   [IMF,RESIDUAL] = VMD(X) returns intrinsic mode functions (IMFs) and a
%   residual signal corresponding to the variational mode decomposition
%   (VMD) of X, with default decomposition parameters. X can be a vector or
%   a timetable with a single variable containing a vector. X must be a
%   real signal. When X is a vector, IMF is a matrix where each column
%   stores an extracted intrinsic mode function, and RESIDUAL is a column
%   vector storing the residual. When X is a timetable, IMF is a timetable
%   with multiple single variables where each variable stores a mode
%   function, and RESIDUAL is a timetable with a single variable. X can be
%   double or single precision.
%
%   [IMF,RESIDUAL] = VMD(X,'Name1',Value1,'Name2',Value2,...)
%   specifies name-value pairs that configure the initial optimization
%   settings and the decomposition stopping criterion to be used for VMD.
%   The supported name-value pairs are:
%
%      'NumIMFs':            Number of decomposition IMFs. The default
%                            value is 5.
%
%      'MaxIterations':      Maximum number of optimization iterations. The
%                            optimization process stops when the number of
%                            iterations is greater than MaxIterations. The
%                            default value is 500.
%
%      'AbsoluteTolerance':  Mode convergence tolerance. The optimization
%      'RelativeTolerance':  process stops when two conditions are met at
%                            the same time: 1) The average squared absolute
%                            improvement toward convergence of IMFs in two
%                            consecutive iterations is less than
%                            AbsoluteTolerance; and 2) The average relative
%                            improvement toward convergence of IMFs in two
%                            consecutive iterations is less than
%                            RelativeTolerance. AbsoluteTolerance and
%                            RelativeTolerance are both specified as
%                            positive real values and default to 5e-6 and
%                            AbsoluteTolerance*1e3, respectively.
% 
%      'CentralFrequencies': Initial central IMF frequencies, specified as
%                            a vector of length NumIMFs. Vector values must
%                            be within the range, [0,0.5] cycles/sample,
%                            which is equivalent to the frequency range,
%                            [0,pi] radians/sample.
%
%      'InitializeMethod':   Method to initialize the central frequencies.
%                            'CentralFrequencies' and 'InitializeMethod'
%                            are mutually exclusive. The methods can be:
%
%                            'random' - Initialize the central frequencies
%                            as random numbers distributed uniformly in the
%                            interval [0,0.5].
%
%                            'grid' - Initialize the central frequencies as
%                            a uniformly sampled grid in the interval
%                            [0,0.5].
%                            
%                            'peaks' - Initialize the central frequencies
%                            as the peak locations of the signal in the
%                            frequency domain (default).
%
%      'InitialIMFs':        Initial IMFs, specified as a real matrix with
%                            rows corresponding to time samples and columns
%                            corresponding to modes. The default value is a
%                            matrix of zeros.
%
%      'PenaltyFactor':      Penalty factor for reconstruction fidelity,
%                            specified as a positive real value. The
%                            smaller its the value, the stricter the data
%                            fidelity. The default value is 1000.
%
%      'InitialLM':          Initial frequency-domain Lagrange multiplier
%                            over the interval, [0,0.5]. The multiplier
%                            enforces the reconstruction constraint. The
%                            default value is a complex vector of zeros.
%                            The length of the multiplier depends on the
%                            input size. See documentation for more
%                            details.
%
%      'LMUpdateRate':       Update rate for the Lagrange multiplier in
%                            each iteration. A higher rate results in
%                            faster convergence but increases the chance of
%                            getting stuck in a local optimum. The default
%                            value is 0.01.
%
%      'Display':            Set to true to display the average absolute
%                            and relative improvement of modes and central
%                            frequencies every 20 iterations, and show the
%                            final stopping information. The default is
%                            false.
%
%   [IMF,RESIDUAL,INFO] = VMD(...) returns the IMFs, the residual, and a
%   structure containing these fields:
%      
%      ExitFlag:            Termination flag. A value of 0 indicates the
%                           algorithm stopped when it reached the maximum
%                           number of iterations. A value of 1 indicates
%                           the algorithm stopped when it met the absolute
%                           and relative tolerances.
%
%      CentralFrequencies:  Central frequencies of the IMFs.
%
%      NumIterations:       Total number of iterations.
%
%      AbsoluteImprovement: Average squared absolute improvement toward
%                           convergence of the IMFs between the final two
%                           iterations.
%
%      RelativeImprovement: Average relative improvement toward
%                           convergence of the IMFs between the final two
%                           iterations.
%
%      LagrangeMultiplier:  Frequency-domain Lagrange multiplier at the 
%                           last iteration.
%
%   VMD(...) with no output arguments plots the original signal, the 
%   residual signal, and the IMFs in the same figure.
%
%   % EXAMPLE 1:
%      % Compute and display the VMD of a signal
%      t = 0:1e-3:1;
%      x1 = cos(2*pi*2*t);
%      x2 = 1/4*cos(2*pi*24*t);
%      x3 = 1/16*cos(2*pi*288*t);
%      x = x1 + x2 + x3 + 0.1*randn(1,length(t));
%      vmd(x,'NumIMFs',3,'Display',true)
%
%   % EXAMPLE 2:
%      % Compute the VMD of a signal and output decomposition details. 
%      % Check that the summation of the IMFs and the residual returns the
%      % original signal
%      t = 0:1e-3:4;
%      x1 = sin(2*pi*50*t) + sin(2*pi*200*t);
%      x2 = sin(2*pi*25*t) + sin(2*pi*100*t) + sin(2*pi*250*t);
%      x = [x1 x2] + 0.1*randn(1,length(t)*2);
%      [IMFs,residual,info] = vmd(x,'MaxIterations',600);
%      max(x(:) - (sum(IMFs,2)+residual))
%
%   See also HHT and EMD.
%

%   Copyright 2019 The MathWorks, Inc.

%#codegen

%---------------------------------
% Check inputs/outputs
narginchk(1,23);
if coder.target('MATLAB') % for MATLAB
    nargoutchk(0,3);
else
    nargoutchk(1,3);
end

% Parse input
[x,td,isTT] = parseInput(x);

% Parse name-value pairs
opts = signalwavelet.internal.vmd.vmdParser(length(x),class(x),varargin{:});

[IMF,residual,info] = computeVMD(x,opts);

if (nargout == 0) && coder.target('MATLAB')
    signalwavelet.internal.convenienceplot.imfPlot(x,IMF,residual,td,'vmd');
end

if isTT && coder.target('MATLAB')
    IMF = array2timetable(IMF,'RowTimes',td);
    residual = array2timetable(residual,'RowTimes',td);
end

if nargout > 0
    varargout{1} = IMF;
end

if nargout > 1
    varargout{2} = residual;
end

if nargout > 2
    varargout{3} = info;
end           
end

function [IMFs,residual,info] = computeVMD(data,opts)
% VMD process
x = cast(data,opts.DataType);

nfft = opts.FFTLength;
penaltyFactor = opts.PenaltyFactor;
numIMFs = opts.NumIMFs;
relativeDiff = cast(inf,opts.DataType);
absoluteDiff = relativeDiff;
tau = opts.LMUpdateRate; % Lagrange multiplier update rate

% Reduce edge effect by mirroring signal; 0 stands for apply mirror
% signal frequency domain with full bandwidth
sigFDFull = signalBoundary(x,opts,0); 
% Get half of the bandwidth
sigFD = sigFDFull(1:opts.NumHalfFreqSamples);

% fft for initial IMFs and get half of bandwidth
initIMFfdFull = fft(opts.InitialIMFs,nfft);
initIMFfd = initIMFfdFull(1:opts.NumHalfFreqSamples,:) + eps; 
IMFfd = initIMFfd;
sumIMF = sum(IMFfd,2);
LM = opts.InitialLM(:); % Lagrange Multiplier

% Frequency vector from [0,0.5) for odd nfft and [0,0.5] for even nfft
f = cast(((0:(nfft/2))/nfft).',opts.DataType);
% Get the initial central frequencies
if strcmp(opts.InitializeMethod,'peaks')
    centralFreq = initialCentralFreqByFindPeaks(abs(sigFD),f,opts);
else
    centralFreq = opts.CentralFrequencies(:);
end

% Progress display set-up
if coder.target('MATLAB') && opts.Display
    fprintf('#Iteration  |  Absolute Improvement  |  Relative Improvement  |  Central Frequencies \n');
    formatstr = '  %5.0f     |       %8.4e       |        %8.4e      |  %s \n';
end
iter = 0;

initIMFNorm = abs(initIMFfd).^2;
normIMF = zeros(size(initIMFfd,1),size(initIMFfd,2),opts.DataType);
%% Optimization iterations    
while (iter < opts.MaxIterations && (relativeDiff > opts.RelativeTolerance ||...
        absoluteDiff > opts.AbsoluteTolerance))   
    for kk = 1:numIMFs
       sumIMF = sumIMF - IMFfd(:,kk);          
       IMFfd(:,kk) = (sigFD - sumIMF + LM/2)./...
           (1+penaltyFactor*(f - centralFreq(kk)).^2);
       normIMF(:,kk) = abs(IMFfd(:,kk)).^2;
       centralFreq(kk) = (f.'*normIMF(:,kk))/sum(normIMF(:,kk));           
       sumIMF = sumIMF + IMFfd(:,kk);            
    end

    LM = LM + tau*(sigFD-sumIMF);
    absDiff = mean(abs(IMFfd-initIMFfd).^2);
    absoluteDiff = sum(absDiff);
    relativeDiff = sum(absDiff./mean(initIMFNorm));

    % Sort IMF and central frequecies in descend order
    % In ADMM, the IMF with greater power will be substracted first
    [~,sortedIndex] = sort(sum(abs(IMFfd).^2),'descend');
    IMFfd = IMFfd(:,sortedIndex);
    centralFreq = centralFreq(sortedIndex(1:length(centralFreq)));
    initIMFfd = IMFfd;
    initIMFNorm = normIMF;  
    iter = iter + 1;
    
    % Progress display
    if coder.target('MATLAB') && opts.Display && (iter ==1 || ~mod(iter,20)...
            || iter >= opts.MaxIterations || (relativeDiff <= opts.RelativeTolerance...
            && absoluteDiff <= opts.AbsoluteTolerance))
        fprintf(formatstr, iter, absoluteDiff, relativeDiff, mat2str(centralFreq.',3));
    end
end
%% Convert to time domain signal     
% Transform to time domain
IMFfdFull = cast(complex(zeros(nfft,numIMFs)),opts.DataType);
IMFfdFull(1:size(IMFfd,1),:) = IMFfd;
if ~mod(opts.FFTLength,2)
    IMFfdFull(size(IMFfd,1)+1:end,:) = conj(IMFfd(end-1:-1:2,:));
else
    IMFfdFull(size(IMFfd,1)+1:end,:) = conj(IMFfd(end:-1:2,:));
end

[~,index] = sort(centralFreq,'descend');   
IMFs = signalBoundary(IMFfdFull(:,index),opts,1); % convert modes back to time domain

%% Output information
info = struct('ExitFlag',0,...
    'CentralFrequencies',centralFreq(index),...
    'NumIterations',iter,...
    'AbsoluteImprovement',absoluteDiff,...
    'RelativeImprovement',relativeDiff,...
    'LagrangeMultiplier',complex(LM));
% Specify stopping flag
if iter < opts.MaxIterations
    info.ExitFlag = 1;
end

% Calculate residual
residual = x - sum(IMFs,2);

% Show stopping information
if coder.target('MATLAB') && opts.Display
    switch info.ExitFlag
        case 0
            finalDisplayStr = getString(message('shared_signalwavelet:vmd:vmd:MaxNumIterationHit', opts.MaxIterations));
        case 1
            finalDisplayStr = getString(message('shared_signalwavelet:vmd:vmd:ToleranceHit', num2str(opts.AbsoluteTolerance),num2str(opts.RelativeTolerance)));
    end
    disp([newline finalDisplayStr newline]);
end
end

function y = signalBoundary(x,opts,isInverse)
% signalBoundary applies mirroring to signal if ifInverse is 0 and removes
% mirrored signal otherwise. Mirror extension of the signal by half its
% length on each side. Removing mirrored signal is a inverse process of the
% mirror extension.
if isInverse % removed mirrored signal
    xr = real(ifft(x,opts.FFTLength));
    y = xr(opts.HalfSignalLength+1:opts.MirroredSignalLength-opts.HalfSignalLength,:);
else % apply mirror to signal
    xr = [x(opts.HalfSignalLength:-1:1); x;...
        x(opts.SignalLength:-1:ceil(opts.SignalLength/2)+1)];
    y = fft(xr,opts.FFTLength);
end
end

function centralFreq = initialCentralFreqByFindPeaks(x,f,opts)  
% Initialize central frequencies by finding the locations of signal peaks
% in frequency domain by using findpeaks function. The number of peaks is
% determined by NumIMFs.
BW = cast(2/opts.FFTLength,opts.DataType); % bandwidth of signal
minBWGapIndex = 2*BW/f(2);

x(x<mean(x)) = mean(x);
TF = islocalmax(x,'MinSeparation',minBWGapIndex);
pkst = x(TF);
locst = f(TF);
numpPeaks = length(pkst);

% Check for DC component
if x(1) >= x(2)
    pks = zeros(numpPeaks+1,1);
    locs = pks;
    pks(2:length(pkst)+1) = pkst;
    locs(2:length(pkst)+1) = locst;
    pks(1) = x(1);
    locs(1) = f(1);
else
    pks = zeros(numpPeaks,1);
    locs = pks;
    pks(1:length(pkst)) = pkst;
    locs(1:length(pkst)) = locst;
end   

[~,index] = sort(pks,'descend');
centralFreq = 0.5*rand(opts.NumIMFs,1,opts.DataType);

% Check if the number of peaks is less than number of IMFs
if length(locs) < opts.NumIMFs
    centralFreq(1:length(locs(index))) = locs;
else
    centralFreq(1:opts.NumIMFs) = locs(index(1:opts.NumIMFs));
end
end

function [x,td,isTT] = parseInput(x)
% Parse the input    
isTT = isa(x,'timetable');
if isTT
    % Verify class if not timetable in codegen
    coder.internal.errorIf(~coder.target('MATLAB'),'shared_signalwavelet:vmd:vmd:TimetableNotSupportedCodegen');
        
    signalwavelet.internal.util.utilValidateattributesTimetable(x,...
        {'sorted','singlechannel','regular'});
    [x, ~, td] = signalwavelet.internal.util.utilParseTimetable(x);   
    validateattributes(x, {'single','double',},{'nonnan',...
        'finite','real', 'nonsparse','vector'},'vmd','the variable of the timetable X');
else
    td = cast((1:length(x)).',class(x));
    % Keep the timetable check so that the error message lists all valid
    % data types.
    validateattributes(x, {'single','double','timetable'},{'nonnan',...
            'finite','real', 'nonsparse','vector'},'vmd','X');
end
coder.internal.errorIf(length(x)<2,'shared_signalwavelet:vmd:vmd:InvalidNumDataSamples');
x = x(:);
end
