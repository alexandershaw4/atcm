function [Y,w,Fspec,units,MAG,PHA] = Alex_FourierTF(P,M,U,opts)
% Alex_FourierTF
% Straight Fourier / simulation-based spectral drop-in for DCM/TCM models.
%
% This is deliberately different from Alex_LaplaceTFwDNew:
%   - no resolvent inv(sI - A) is used
%   - no explicit eigenmode or Koopman/DMD decomposition is used
%   - the nonlinear model is simulated in the time domain
%   - observed regional signals are Fourier transformed and interpolated to M.Hz
%
% Intended use:
%   [Y,w,Fspec,units,MAG,PHA] = Alex_FourierTF(P,M,U,opts)
%
% Outputs are shaped to match the Laplace/Eigenmode/Koopman drop-ins:
%   Y      : {CSD}, where CSD is [frequency x source x source]
%   w      : frequencies requested in M.Hz
%   Fspec  : struct containing FFT details and regional spectra
%   units  : simulated states/signals and options
%   MAG    : per-source state FFT magnitudes
%   PHA    : per-source state FFT phases in degrees
%
% Notes:
%   This is a useful diagnostic comparator, but less efficient and usually less
%   smooth than the Laplace-domain transfer function. It captures nonlinear
%   time-domain behaviour, but delay handling depends on whether delays are
%   already present in M.f. The explicit Laplace delay matrix D is not used here.
%
% AS2026

if nargin < 4, opts = struct(); end
if ~isfield(opts,'dt'),          opts.dt = 1e-3; end              % seconds
if ~isfield(opts,'T'),           opts.T = 4.0; end                % seconds
if ~isfield(opts,'burnin'),      opts.burnin = 0.4; end           % seconds
if ~isfield(opts,'nrep'),        opts.nrep = 2; end               % repetitions for averaging
if ~isfield(opts,'x0_scale'),    opts.x0_scale = 1e-4; end
if ~isfield(opts,'noise_scale'), opts.noise_scale = 1e-3; end     % state noise for endogenous drive
if ~isfield(opts,'drive_scale'), opts.drive_scale = 1.0; end
if ~isfield(opts,'demean'),      opts.demean = true; end
if ~isfield(opts,'window'),      opts.window = 'hann'; end
if ~isfield(opts,'eps_reg'),     opts.eps_reg = 1e-12; end
if ~isfield(opts,'use_rk2'),     opts.use_rk2 = true; end
if ~isfield(opts,'seed'),        opts.seed = []; end

if ~isempty(opts.seed)
    rng(opts.seed);
end

if isnumeric(P), P = spm_unvec(P,M.P); end
if isstruct(P) && isfield(P,'p'), P = P.p; end

if isfield(M,'endogenous') && M.endogenous
    Input = 0;
else
    Input = 1;
end

% Optional fixed point, matching the other drop-ins
if isfield(M,'fixedpoint') && M.fixedpoint == 1
    x = atcm.fun.alexfixed(P,M,1e-10,[],[],1000);
    M.x = spm_unvec(x,M.x);
end

w   = M.Hz(:);
x0  = M.x(:);
nx  = numel(x0);
Ns  = size(M.x,1);

% Time base
Ntotal = max(8,round(opts.T / opts.dt));
Nburn  = max(0,round(opts.burnin / opts.dt));
Nkeep  = Ntotal - Nburn;
if Nkeep < 8
    error('Alex_FourierTF:NotEnoughSamples','T - burnin gives too few samples.');
end

Fs = 1 / opts.dt;
f_fft = (0:Nkeep-1)' * Fs / Nkeep;
keep_freq = f_fft <= max(w)*1.25 & f_fft >= 0;
f_use = f_fft(keep_freq);

% Window, avoiding toolbox dependencies
switch lower(opts.window)
    case {'hann','hanning'}
        ww = 0.5 - 0.5*cos(2*pi*(0:Nkeep-1)'/(Nkeep-1));
    case 'hamming'
        ww = 0.54 - 0.46*cos(2*pi*(0:Nkeep-1)'/(Nkeep-1));
    otherwise
        ww = ones(Nkeep,1);
end
ww = ww(:);
win_norm = sqrt(sum(ww.^2));

% Optional external frequency spectrum, used only to shape a synthetic drive
Uomega = ones(numel(w),1);
if isfield(M,'external_spectrum')
    Uomega = M.external_spectrum(:);
end

% Build a simple broadband input in the time domain when exogenous.
% This is intentionally modest: for exact linear FRFs, use the Laplace routine.
make_input = @(tidx) local_get_input(U,tidx,Input,opts.drive_scale);

% Accumulators
PSD_acc = zeros(Ns,numel(w));
CSD_acc = zeros(numel(w),Ns,Ns);
MAG_acc = cell(Ns,1);
PHA_acc = cell(Ns,1);
for ii = 1:Ns
    MAG_acc{ii} = [];
    PHA_acc{ii} = [];
end

all_y = zeros(Nkeep,Ns,opts.nrep);
all_x = zeros(Nkeep,nx,opts.nrep);

for rr = 1:opts.nrep
    x = x0 + opts.x0_scale * randn(size(x0));
    Xkeep = zeros(Nkeep,nx);

    kk = 0;
    for t = 1:Ntotal
        ut = make_input(t);

        dx = feval(M.f, spm_unvec(x,M.x), ut, P, M);
        dx = dx(:);

        % Weak stochastic drive helps reveal endogenous resonances.
        if ~Input && opts.noise_scale > 0
            dx = dx + opts.noise_scale * randn(size(dx));
        end

        if opts.use_rk2
            xmid = x + 0.5 * opts.dt * dx;
            dxmid = feval(M.f, spm_unvec(xmid,M.x), ut, P, M);
            dxmid = dxmid(:);
            if ~Input && opts.noise_scale > 0
                dxmid = dxmid + opts.noise_scale * randn(size(dxmid));
            end
            x = x + opts.dt * dxmid;
        else
            x = x + opts.dt * dx;
        end

        if any(~isfinite(x))
            warning('Alex_FourierTF:NonFiniteState','Non-finite state encountered; clipping and continuing.');
            x(~isfinite(x)) = 0;
            x = max(min(x,1e6),-1e6);
        end

        if t > Nburn
            kk = kk + 1;
            Xkeep(kk,:) = x(:).';
        end
    end

    if opts.demean
        Xkeep = Xkeep - mean(Xkeep,1);
    end

    Yreg = zeros(Nkeep,Ns);
    for ii = 1:Ns
        win = ii:Ns:nx;
        Cw  = exp(P.J(win));

        % Dynamic observer is not used here because it is frequency-domain.
        % Use the static observer weights to create a time-domain regional signal.
        Yreg(:,ii) = Xkeep(:,win) * Cw(:);

        if isfield(P,'L') && numel(P.L) >= ii
            Yreg(:,ii) = exp(P.L(ii)) * Yreg(:,ii);
        end
    end

    if opts.demean
        Yreg = Yreg - mean(Yreg,1);
    end

    all_y(:,:,rr) = Yreg;
    all_x(:,:,rr) = Xkeep;

    % Regional FFT and cross-spectra
    Yfft = fft(Yreg .* ww,[],1) / max(win_norm,opts.eps_reg);
    Yfft = Yfft(keep_freq,:);

    for ii = 1:Ns
        PSD_i = abs(Yfft(:,ii)).^2;
        PSD_acc(ii,:) = PSD_acc(ii,:) + interp1(f_use,PSD_i,w,'linear','extrap').';
    end

    for ii = 1:Ns
        for jj = 1:Ns
            Cij = Yfft(:,ii) .* conj(Yfft(:,jj));
            CSD_acc(:,ii,jj) = CSD_acc(:,ii,jj) + interp1(f_use,Cij,w,'linear','extrap');
        end
    end
    
    % Smooth magnitudes (keeps behaviour of |.| then smooth)
    dw = mean(diff(w));
    if Ns == 1
        CSD = atcm.fun.agauss_smooth(abs(CSD_acc), dw * exp(P.d(3)));
    else
        for ii = 1:Ns
            for jj = 1:Ns
                CSD_acc(:,ii,jj) = atcm.fun.agauss_smooth(abs(CSD_acc(:,ii,jj)), dw * exp(P.d(3)));
                CSD(:,jj,ii) = CSD(:,ii,jj);
            end
        end
    end


    % State-level FFTs per source for MAG/PHA compatibility
    Xfft = fft(Xkeep .* ww,[],1) / max(win_norm,opts.eps_reg);
    Xfft = Xfft(keep_freq,:);
    for ii = 1:Ns
        win = ii:Ns:nx;
        Xfi = zeros(numel(w),numel(win));
        for kk2 = 1:numel(win)
            Xfi(:,kk2) = interp1(f_use,Xfft(:,win(kk2)),w,'linear','extrap');
        end
        if rr == 1
            MAG_acc{ii} = abs(Xfi.');
            PHA_acc{ii} = angle(Xfi.') * 180/pi;
        else
            MAG_acc{ii} = MAG_acc{ii} + abs(Xfi.');
            PHA_acc{ii} = PHA_acc{ii} + angle(Xfi.') * 180/pi;
        end
    end
end

PSD = PSD_acc / opts.nrep;
CSD = CSD_acc / opts.nrep;

MAG = cell(Ns,1);
PHA = cell(Ns,1);
for ii = 1:Ns
    MAG{ii} = MAG_acc{ii} / opts.nrep;
    PHA{ii} = PHA_acc{ii} / opts.nrep;
end

% Match your existing convention: smooth magnitudes and return real CSD.
dw = mean(diff(w));
if isfield(P,'d') && numel(P.d) >= 3
    if Ns == 1
        CSD = atcm.fun.agauss_smooth(abs(CSD), dw * exp(P.d(3)));
    else
        for ii = 1:Ns
            for jj = 1:Ns
                CSD(:,ii,jj) = atcm.fun.agauss_smooth(abs(CSD(:,ii,jj)), dw * exp(P.d(3)));
                CSD(:,jj,ii) = CSD(:,ii,jj);
            end
        end
    end
else
    CSD = abs(CSD);
end

Y = {CSD};

Fspec = struct();
Fspec.f_fft = f_use;
Fspec.PSD = PSD;
Fspec.Input = Input;
Fspec.Uomega = Uomega;
Fspec.window = opts.window;

units = struct();
units.x0 = x0;
units.freq = w;
units.dt = opts.dt;
units.Fs = Fs;
units.T = opts.T;
units.burnin = opts.burnin;
units.nrep = opts.nrep;
units.t = (0:Nkeep-1)' * opts.dt;
units.X = all_x;
units.LFP = all_y;
units.opts = opts;

end

function ut = local_get_input(U,tidx,Input,drive_scale)
if ~Input
    ut = 0;
    return;
end

ut = drive_scale;
if nargin < 1 || isempty(U)
    return;
end

try
    if isstruct(U) && isfield(U,'u') && ~isempty(U.u)
        if size(U.u,2) >= tidx
            ut = U.u(:,tidx) * drive_scale;
        else
            ut = U.u(:,end) * drive_scale;
        end
    elseif isnumeric(U) && ~isempty(U)
        if size(U,2) >= tidx
            ut = U(:,tidx) * drive_scale;
        else
            ut = U(:,end) * drive_scale;
        end
    end
catch
    ut = drive_scale;
end
end
