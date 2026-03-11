function [Y,w,Kmodes,units,MAG,PHA] = Alex_KoopmanTF(P,M,U,opts)
% Robust Koopman / DMD spectral surrogate for nonlinear DCM/TCM
%
% Uses:
%   - simulation of nonlinear dynamics
%   - SVD-truncated DMD
%   - eigenvalue-wise continuous-time conversion
%   - modal spectral reconstruction
%
% Optional opts fields:
%   .dt        simulation step (default 1e-3)
%   .T         number of simulation samples (default 4000)
%   .rank      DMD truncation rank
%   .demean    demean snapshots (default true)
%   .eps_reg   small regulariser (default 1e-8)
%   .x0_scale  perturbation scale around initial condition (default 1e-4)

if nargin < 4, opts = struct(); end
if ~isfield(opts,'dt'),       opts.dt = 1e-3; end
if ~isfield(opts,'T'),        opts.T = 4000; end
if ~isfield(opts,'demean'),   opts.demean = true; end
if ~isfield(opts,'eps_reg'),  opts.eps_reg = 1e-8; end
if ~isfield(opts,'x0_scale'), opts.x0_scale = 1e-4; end

if isnumeric(P), P = spm_unvec(P,M.P); end
if isstruct(P) && isfield(P,'p'), P = P.p; end

if isfield(M,'fixedpoint') && M.fixedpoint == 1
    x = atcm.fun.alexfixed(P,M,1e-10,[],[],1000);
    M.x = spm_unvec(x,M.x);
end

w  = M.Hz(:);
x0 = M.x(:);
n  = numel(x0);
Ns = size(M.x,1);

% Slight perturbation so snapshots are not perfectly rank-deficient
x = x0 + opts.x0_scale * randn(size(x0));

X  = zeros(n, opts.T-1);
Xp = zeros(n, opts.T-1);

for t = 1:(opts.T-1)
    ut = 0;
    if nargin >= 3 && ~isempty(U)
        try
            if isstruct(U) && isfield(U,'u') && numel(U.u) >= t
                ut = U.u(:,t);
            elseif isnumeric(U) && numel(U) >= t
                ut = U(:,min(t,size(U,2)));
            end
        catch
            ut = 0;
        end
    end

    dx = feval(M.f, spm_unvec(x,M.x), ut, P, M);
    dx = dx(:);

    % RK2 / midpoint step: more stable than Euler
    xmid = x + 0.5 * opts.dt * dx;
    dxmid = feval(M.f, spm_unvec(xmid,M.x), ut, P, M);
    dxmid = dxmid(:);

    xnext = x + opts.dt * dxmid;

    X(:,t)  = x;
    Xp(:,t) = xnext;
    x       = xnext;
end

% optional demeaning
xmean = zeros(n,1);
if opts.demean
    xmean = mean(X,2);
    X  = X  - xmean;
    Xp = Xp - mean(Xp,2);
end

% SVD truncation
[Uu,S,Vv] = svd(X,'econ');
sing = diag(S);

if isfield(opts,'rank') && ~isempty(opts.rank)
    r = min(opts.rank, numel(sing));
else
    energy = cumsum(sing.^2) / max(sum(sing.^2), eps);
    r = find(energy >= 0.999, 1, 'first');
    if isempty(r), r = min([size(Uu,2), 12]); end
end

Ur = Uu(:,1:r);
Sr = S(1:r,1:r);
Vr = Vv(:,1:r);

Atilde = Ur' * Xp * Vr / (Sr + opts.eps_reg*eye(r));
[W,L]  = eig(Atilde);

lambda_d = diag(L);
mu_c     = log(lambda_d) / opts.dt;   % continuous-time poles, eigenvalue-wise
Phi      = Xp * Vr / (Sr + opts.eps_reg*eye(r)) * W;  % DMD modes

% modal amplitudes from initial condition
b = pinv(Phi) * (x0 - xmean);

PSD = zeros(Ns, numel(w));
MAG = cell(Ns,1);
PHA = cell(Ns,1);

Kmodes = struct();
Kmodes.lambda_d = lambda_d;
Kmodes.mu_c     = mu_c;
Kmodes.Phi      = Phi;
Kmodes.b        = b;
Kmodes.rank     = r;
Kmodes.Atilde   = Atilde;
Kmodes.singular_values = sing;

for ii = 1:Ns
    win = ii:Ns:n;
    Cw  = exp(P.J(win));
    MG  = complex(zeros(numel(win), numel(w)));
    y   = complex(zeros(1, numel(w)));

    % observation weights on each mode
    alpha = zeros(r,1);
    for k = 1:r
        alpha(k) = (Cw(:).' * Phi(win,k)) * b(k);
    end

    for j = 1:numel(w)
        s = 1i * 2*pi * w(j);

        Hj = 0;
        xj = zeros(n,1);
        for k = 1:r
            denom = (s - mu_c(k));
            if abs(denom) < 1e-8
                denom = denom + 1e-8;
            end
            coeff = b(k) / denom;
            xj = xj + Phi(:,k) * coeff;
            Hj = Hj + alpha(k) / denom;
        end

        MG(:,j) = xj(win);
        y(j)    = Hj;
    end

    if isfield(P,'L') && numel(P.L) >= ii
        PSD(ii,:) = exp(P.L(ii)) * y;
    else
        PSD(ii,:) = y;
    end

    MAG{ii} = abs(MG);
    PHA{ii} = angle(MG) * 180/pi;
end

CSD = zeros(numel(w), Ns, Ns);
for i = 1:Ns
    CSD(:,i,i) = PSD(i,:).';
    for j = 1:Ns
        if i ~= j
            CSD(:,i,j) = PSD(i,:).' .* conj(PSD(j,:).');
            CSD(:,j,i) = CSD(:,i,j);
        end
    end
end

% optional smoothing in your style
if isfield(P,'d') && numel(P.d) >= 3
    dw = mean(diff(w));
    for i = 1:Ns
        for j = 1:Ns
            CSD(:,i,j) = atcm.fun.agauss_smooth(abs(CSD(:,i,j)), dw * exp(P.d(3)));
        end
    end
else
    CSD = abs(CSD);
end

Y = {CSD};

units = struct();
units.X = X;
units.Xp = Xp;
units.xmean = xmean;
units.freq = w;
units.rank = r;
units.x0 = x0;