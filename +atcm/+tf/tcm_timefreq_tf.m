function [TFR, t, f, Y, WT] = tcm_timefreq_tf(P, M, U, tfopts)
% TCM_TIMEFREQ_TF  Time–frequency forward operator for a linearised TCM.
%
%   [TFR,t,f,Y,WT] = tcm_timefreq_tf(P,M,U,tfopts)
%
% Steps:
%   1) Linearise TCM: [A,B,C,Delays] = atcm.linearise(P,M)
%   2) Generate time-domain output y(t) to input U(t) with delays
%   3) Compute complex Morlet TFR: WT = cwt(y, 'amor', fs)
%   4) Baseline normalise to % change: TFR (real-valued)
%
% Inputs
%   P      : parameter struct
%   M      : model struct; must contain fs, T, output index, etc.
%   U      : struct with U.t (s × 1), U.u (nu × s) input
%   tfopts : struct with fields:
%              .fband = [4 120];           % frequency range (Hz)
%              .nb = 48;                   % # wavelet voices (or use default)
%              .baseline = [-0.5 0];       % baseline window (s)
%              .output_ix = 1;             % which sensor/virtual channel
%              .dmethod = 'shift';         % delays: 'shift' or 'pade'
%
% Outputs
%   TFR : (nf × nt) % change
%   t   : time vector (s)
%   f   : frequency vector (Hz)
%   Y   : time-domain channel output
%   WT  : complex CWT coefficients (nf × nt)

arguments
    P
    M
    U
    tfopts;
    %tfopts.fband (1,2) double = [4 120]
    %tfopts.nb (1,1) double = 0         % 0 -> let cwt choose
    %tfopts.baseline (1,2) double = [-0.5 0]
    %tfopts.output_ix (1,1) double = 1
    %tfopts.dmethod (1,1) string = "shift"
end

if nargin < 4 || isempty(tfopts); tfopts = struct; end

if ~isfield(tfopts,'fband'); tfopts.fband = [1 90]; end
if ~isfield(tfopts,'nb');    tfopts.nb = 0; end
if ~isfield(tfopts,'baseline'); tfopts.baseline = [-1000 1]; end
if ~isfield(tfopts,'output_ix'); tfopts.output_ix = 1; end
if ~isfield(tfopts,'dmethod'); tfopts.dmethod = 'shift'; end



fs  = M.fs;                 % sampling rate (Hz)
%t   = 0:1/fs:M.T(end);           % simulation horizon (s)
t = M.T;
nt  = numel(t);

% ----- 1) Linearise the TCM (your existing code) -------------------------
% Expect: A (nx×nx), B (nx×nu), C (ny×nx), delays struct with .Ax/.Cx/.Ux
[A,B,C,Del] = atcm.linearise(P,M,U);    %#ok<NASGU>  % you already have this

% Ensure B has at least one input column
if isempty(B)
    B = zeros(size(A,1),1);
end
nu = size(B,2);

% --- build input over t (self-consistent with nu) ---
if ~(exist('U','var') && isstruct(U) && isfield(U,'t') && isfield(U,'u') && ~isempty(U.u))
    uu = zeros(nu, nt);
    uu(:, t>=0) = 1;                      % step input at 0 s
else
    uu = interp1(U.t(:), U.u.', t, 'previous', 'extrap').';
    % match input dimensionality expected by B
    if size(uu,1) ~= nu
        uu = uu(1:min(end,nu),:);
        if size(uu,1) < nu
            uu = [uu; zeros(nu - size(uu,1), nt)];
        end
    end
end

% --- delay sanitiser (robust to NaN, vectors, negatives) ---
fs_ = M.fs;
dU = 0; dC = 0;

if exist('Del','var') && isstruct(Del)
    % U delay
    if isfield(Del,'U') && ~isempty(Del.U)
        dUval = Del.U;
        if ~isscalar(dUval), dUval = max(dUval(:)); end
        if isfinite(dUval) && dUval > 0
            dU = round(dUval * fs_);
        end
    end
    % C delay
    if isfield(Del,'C') && ~isempty(Del.C)
        dCval = Del.C;
        if ~isscalar(dCval), dCval = max(dCval(:)); end
        if isfinite(dCval) && dCval > 0
            dC = round(dCval * fs_);
        end
    end
end

% Final guards
dU = max(0, floor(dU));
dC = max(0, floor(dC));

% --- simulate ---
nx = size(A,1);
x  = zeros(nx,1);
Y  = zeros(1,nt);

for k = 2:nt
    dt = 1/fs_;
    ku = max(1, k - dU);                  % always >=1 integer
    x  = x + dt*(A*x + B*uu(:,ku));
    yk = C'*x;                             % handle multi-row C if needed
    Y(k) = yk(1);
end

% Optionally subtract mean pre-stim
if any(t<0)
    Y = Y - mean(Y(t<0));
end

% ----- 3) Complex Morlet CWT ---------------------------------------------
if tfopts.nb>0
    [WT,f] = cwt(Y, 'amor', fs, 'FrequencyLimits', tfopts.fband, 'VoicesPerOctave', tfopts.nb);
else
    [WT,f] = cwt(Y, 'amor', fs, 'FrequencyLimits', tfopts.fband);
end
% WT is nf × nt complex; power:
Pwr = abs(WT).^2;

% ----- 4) Baseline → % change -------------------------------------------
bidx = t>=tfopts.baseline(1) & t<=tfopts.baseline(2);
Pb   = mean(Pwr(:,bidx), 2) + eps;
TFR  = 100*(Pwr./Pb - 1);
end
