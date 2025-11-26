function [TFR, f, t, TFR_trials] = dcm_data_tfm(DCM, fband, baseline)
% DCM → Time–Frequency Map (data)
%
% Inputs
%   DCM      : your DCM struct with xY fields
%   fband    : [fmin fmax] Hz (default: [min(Hz) max(Hz)])
%   baseline : [t0 t1] ms window for baseline (default: min(pst) to 0)
%
% Outputs
%   TFR        : (nf × nt) % change, trial-averaged
%   f          : frequency vector (Hz)
%   t          : time vector (s)
%   TFR_trials : (nf × nt × ntr) % change per trial (optional)
%
%  [TFR, f, t, TFR_trials] = atcm.tf.dcm_data_tfm(DCM,[10 80],[-100 0]);
%
% AS

xY = DCM.xY;
Ytr = xY.series{1};             % trials x samples
fs  = 1./xY.dt;
t   = xY.pst(:)'/1000;          % seconds
if nargin < 2 || isempty(fband), fband = [min(xY.Hz) max(xY.Hz)]; end
if nargin < 3 || isempty(baseline), baseline = [min(xY.pst) 0]; end

% Use MATLAB's cwt with analytic Morlet ('amor') on the first pass.
% We'll match your target band and then select those freqs.
[T1,f] = cwt(Ytr(1,:), 'amor', fs, 'FrequencyLimits', fband);
nf = numel(f); nt = size(T1,2); ntr = size(Ytr,1);

WT = complex(zeros(nf, nt, ntr));
WT(:,:,1) = T1;
for tr = 2:ntr
    WT(:,:,tr) = cwt(Ytr(tr,:), 'amor', fs, 'FrequencyLimits', fband);
end

Pwr = abs(WT).^2;                               % power
bmsk = t >= baseline(1)/1000 & t <= baseline(2)/1000;
Pb   = mean(Pwr(:,bmsk,:), 2) + eps;            % nf × 1 × ntr
TFR_trials = 100*(Pwr./Pb - 1);                 % % change per trial

TFR = mean(TFR_trials, 3);                      % trial average
end
