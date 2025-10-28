function [Y,w,G,units,MAG,PHA] = Alex_LaplaceTFwDNew(P,M,U)
% linearisation and numerical (Laplace) transfer function for a DCM
% witten by Alex Shaw; this version has proper (Laplace domain) handling of
% delays...
%
% takes a dynamical model and observation function;
%
%   dx = f(x,u,P,M)
%    y = g(x,P)
%
% and linearises flow of states [A], inputs [B] and observation [C];
%
%   dx = Ax  + Bu;
%    y = Cx [+ Du];
%
% having computed the linearisation matrices, computes the frequency response 
% at frequencies of interest (vector stored in M.Hz) using the Laplace transform.
%
%   Y(s) = C*inv(sI - A)*BU or  C*inv(sI - A)*x0
%
% Usage: [Y,w,G,units] = atcm.fun.Alex_LaplaceTFwD(P,M,U); 
%
% add sub-structure M.sim with M.sim.pst and M.sim.dt to force the routine
% to also recronstruct the time-domain simulation and return in the 4th
% output, since we can access the magnitude and phase data of each component
% from the Laplace transform.
%
% - Delays enter as elementwise factors: A_eff(s) = A ∘ exp(-s*D)
% - Exogenous case: y(ω) = C * (sI - A_eff)^(-1) * B * u(ω)
% - Endogenous case: S_y(ω) = H(ω) Σ_w H(ω)^* , H(ω) = C * (sI - A_eff)^(-1)
%
% Notes:
%   * If present, P.q (per-region log process noise) sets Σ_w = exp(P.q(i)) I.
%   * P.d(2) can act as a small real shift of s (damping/jitter). Set to 0 if unused.
%
% AS2023 updated Oct 2025

if isnumeric(P), P = spm_unvec(P,M.P); end
if isstruct(P) && isfield(P,'p'), P = P.p; end

if isfield(M,'endogenous') && M.endogenous
    Input = 0;
else
    Input = 1;
end

% Fixed point (optional)
if isfield(M,'fixedpoint') && M.fixedpoint == 1
    x = atcm.fun.alexfixed(P,M,1e-10,[],[],1000);
    M.x = spm_unvec(x,M.x);
end

w   = M.Hz(:);
x0  = M.x(:);                 % not used as drive any more
Ns  = size(M.x,1);

% Linearisation
[~,A,D] = feval(M.f,M.x,0,P,M);         % A, D (states×states)
A       = denan(A);
Bfull   = spm_diff(M.f,M.x,1,P,M,2);    % input Jacobian wrt u
Bfull   = denan(Bfull);

% Default drive spectrum (flat)
Uomega = ones(numel(w),1);
if isfield(M,'external_spectrum')
    Uomega = M.external_spectrum(:);
end

% Optional tiny real shift in s (numerical damping)
damp = 0;
if isfield(P,'d') && numel(P.d) >= 2
    damp = exp(P.d(2));
end

% Optional per-region log-noise power (for endogenous)
have_q = isfield(P,'q') && numel(P.q) >= Ns;

% Prealloc
PSD   = zeros(Ns,numel(w));
MAG   = cell(Ns,1);
PHA   = cell(Ns,1);
G     = [];                 % state-space sys object not used here
Hstore = cell(Ns,1);        % for endogenous path: H rows

for ii = 1:Ns
    % states for region ii (block picking every Ns-th state)
    win = ii:Ns:(length(A));
    n   = numel(win);

    AA  = A(win,win);
    BB  = Bfull(win,:);                      % handle multi-input later by weighting
    Cw  = exp(P.J(win));                     % observer weights (vector)
    Cmat = diag(Cw);                         % lfp readout = Cw' * x_win
    X0 = x0(win);

    drive_scale = 1;
    if isfield(P,'C') && numel(P.C) >= ii
        drive_scale = exp(P.C(ii));
    end

    % If multiple inputs, you can choose one column here or a mix:
    % For now assume single effective drive column:
    if size(BB,2) > 1
        % default: sum/average columns (or pick 1st) — adjust as needed
        BB = sum(BB,2);
    end

    % Per-frequency
    MG = complex(zeros(n,numel(w)));
    y  = complex(zeros(1,numel(w)));

    for j = 1:numel(w)
        s   = damp + 1i*2*pi*w(j);
        E   = exp(-1i*2*pi*w(j) * D(win,win));   % elementwise exp(-s*D_ij) with s=jω
        Aef = AA .* E;                           % A_eff(s) = A ∘ e^{-sD}

        % Resolvent
        Jm  = (s*eye(n)) - Aef;

        if Input
            % Exogenous: y(ω) = C * Jm^{-1} * B * u(ω)
            u_j = Uomega(j) * drive_scale;
            Ym  = (Jm \ (BB * u_j)) + (Jm \ X0);
            MG(:,j) = Ym;
            y(j)    = (Cw.' * Ym);
        else
            % Endogenous: store H row; spectra computed after loop
            % H(ω) maps state noise w (n×1) -> output y (scalar):
            % H_row = Cw' * Jm^{-1}    (1×n)
            Hstore{ii}(j,1:n) = (Cw.' / Jm);
            % Keep MG/y placeholders for completeness (not used)
            MG(:,j) = complex(NaN);
            y(j)    = complex(NaN);

            Hrow = (Cw.' / Jm);       % 1×n complex transfer from state-noise -> output
            qpow = (have_q && numel(P.q)>=ii) * exp(P.q(ii)) + (~(have_q && numel(P.q)>=ii)) * 1;

            % Power (scalar): S_y(ω) = q * sum |H|^2
            psd_i = qpow * sum(abs(Hrow).^2, 2);   % (#ω × 1), real

            % If you want something like MG/PHA for plotting contributions:
            MG(:,j)  = Hrow.';                    % n×1 complex gains per state
            PH(:,j) = angle(MG(:,j))*180/pi;             % optional: "amplitude-like" scalar (no true phase)

        end
    end

    if Input
        % Optional taper on output spectrum magnitude
        if isfield(M,'ham') && M.ham
            Hm = hamming(numel(w));
            y  = y(:) .* Hm(:);
        end

        MAG{ii} = MG;
        PHA{ii} = angle(MG)*180/pi;

        % Electrode scaling (log-gain per region)
        Lgain = 1;
        if isfield(P,'L') && numel(P.L) >= ii
            Lgain = exp(P.L(ii));
        end

        PSD(ii,:) = Lgain * (y(:)).';   % complex response per region
    else
        % Endogenous PSD via sandwich: S_y(ω) = H Σ_w H*
        qpow = 1;
        if have_q, qpow = exp(P.q(ii)); end

        Hrw = Hstore{ii};   % (#ω × n)
        % For Σ_w = q * I, the scalar PSD is q * sum |H_row|^2
        psd_i = qpow * sum(abs(Hrw).^2, 2);   % (#ω × 1), real

        % Electrode scaling (as magnitude scaling)
        Lgain = 1;
        if isfield(P,'L') && numel(P.L) >= ii
            Lgain = exp(P.L(ii));
        end
        PSD(ii,:) = (Lgain * psd_i(:)).';     % real, non-negative

        MAG{ii} = MG;
        PHA{ii} = PH;

    end
end

% --- Cross-spectra construction (simple heuristic, as in your original) ---
% For endogenous path, true cross terms need state-noise cross-covariances.
% We keep original heuristic (scaled product of autospectra).
CSD = zeros(numel(w),Ns,Ns);
for i = 1:Ns
    CSD(:,i,i) = PSD(i,:).';  % auto
    for j = 1:Ns
        if i ~= j
            Lc = 1;
            if isfield(P,'Lc') && numel(P.Lc) >= i
                Lc = exp(P.Lc(i));
            end
            % Use product with conjugate for exogenous (complex), or plain product if real
            if Input
                CSD(:,i,j) = Lc * (PSD(i,:).' .* conj(PSD(j,:).'));
            else
                CSD(:,i,j) = Lc * (PSD(i,:).' .* (PSD(j,:).')); % real, phase-free
            end
            CSD(:,j,i) = CSD(:,i,j);
        end
    end
end

% Smooth magnitudes (keeps behaviour of |.| then smooth)
dw = mean(diff(w));
if Ns == 1
    CSD = atcm.fun.agauss_smooth(abs(CSD), dw * exp(P.d(1)));
else

    for i = 1:Ns
        for j = 1:Ns
            CSD(:,i,j) = atcm.fun.agauss_smooth(abs(CSD(:,i,j)), dw * exp(P.d(1)));
            CSD(:,j,i) = CSD(:,i,j);
        end
    end
end

Y     = {(CSD)};
units = [];

% --- Optional time-domain reconstruction (note: ignores delays in time domain) ---
if isfield(M,'sim') && nargout > 3
    pst = M.sim.pst;
    dt  = M.sim.dt;

    % Re-enter to get MAG/PHA (exogenous only)
    if Input
        M2 = rmfield(M,'sim');
        P2 = P; P2.J(P2.J==-1000)=0;
        [~,~,~,~,MAG,PHA] = atcm.fun.Alex_LaplaceTFwD(P2,M2,U);
    end

    if Input
        for k = 1:Ns
            mag = squeeze(MAG{k});
            the = squeeze(PHA{k});
            series = zeros(size(mag,1), size(mag,2), numel(pst));
            for ii_ = 1:size(mag,1)
                for jj_ = 1:size(mag,2)
                    series(ii_,jj_,:) = mag(ii_,jj_) * sin(2*pi*w(jj_)*pst/1000 - the(ii_,jj_));
                end
            end
            S{k}   = squeeze(sum(series,2));
            S{k}   = S{k} + spm_vec(M.x(k,:,:));
            Cw_top = exp(P.J(k:Ns:end));
            LFP(k,:) = (Cw_top)'*S{k};
        end
        units.series = S;
        units.LFP    = LFP;
    else
        units.series = [];
        units.LFP    = [];
    end

    units.dt     = dt;
    units.pst    = pst;
    units.A      = A;
    units.B      = Bfull;
    units.C      = [];   % observer assembled per-region above
    units.D      = D;
    units.xinit  = M.x(:);
    units.freq   = w(:);
end

return;