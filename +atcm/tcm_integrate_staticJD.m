function [X, info] = tcm_integrate_staticJD(f, M, P, x0, t, u, opts)
%TCM_INTEGRATE_STATICJD  Fast TCM/DCM integrator with static J and D (with stabilisers).
%
% [X, info] = tcm_integrate_staticJD(f, M, P, x0, t, u, opts)
%
% New options:
%   .method     : 'rk4' (default) | 'euler' | 'rosenbrock'   % Rosenbrock–Euler (1st order, stiff)
%   .substeps   : 1 (default). For 'rk4'/'euler', split each dt into N micro-steps.
%   .J, .D      : precomputed Jacobian/Delay (n×n). If missing, we call f once with 3 outputs.
%   .u0         : time to sample u for one-shot J,D eval (default t(1))
%   .maxCond    : guard on cond(A) for A = I - D.*J (default 1e8)
%   .autoAlpha  : false (default). If true, scales D -> alpha*D to make max Re(eig(Jeff)) ≤ 0
%   .alpha      : scalar in (0,1]; if provided, scales D (overrides autoAlpha’s result)
%   .maxStep    : [] (default). If set, caps ||Q*f(x,u)|| per micro-step to this value.
%   .finiteCheck: false (default). If true, error when state becomes non-finite.
%
% Outputs (extra fields):
%   info.alpha      : D scaling used
%   info.evJeff     : eigenvalues of Jeff at start (after scaling)
%   info.maxReJeff  : max real part of eig(Jeff)
%   info.maxAbsJeff : max abs eig(Jeff)

    if nargin < 7 || isempty(opts), opts = struct; end
    method      = getdef(opts,'method','rk4');
    maxCond     = getdef(opts,'maxCond',1e8);
    substeps    = getdef(opts,'substeps',1);
    autoAlpha   = getdef(opts,'autoAlpha',false);
    alpha_user  = getdef(opts,'alpha',[]);
    maxStep     = getdef(opts,'maxStep',[]);
    finiteCheck = getdef(opts,'finiteCheck',false);

    t  = t(:).'; T = numel(t);
    if T < 2, error('t must contain at least two time points.'); end
    if any(diff(t) <= 0), error('t must be strictly increasing.'); end
    dt = diff(t);

    if isempty(u)
        t0 = 0.0; dur = 2+P.d(2); fq = 4+P.d(3);   % 20 Hz for 0.5 s
        u  = @(t) (double(t >= t0 & t < t0+dur) .* sin(2*pi*fq*(t - t0)));
    end

    x0 = x0(:);
    n  = numel(x0);

    % ---------- One-shot J,D (or accept user-provided) ----------
    hasJD = isfield(opts,'J') && ~isempty(opts.J) && isfield(opts,'D') && ~isempty(opts.D);
    if hasJD
        J = opts.J; D = opts.D;
    else
        t0 = getdef(opts,'u0', t(1));
        u0 = sample_u(u, t0, t, 1);
        [~, J, D] = f(x0, u0, P, M);
        if isempty(J), J = zeros(n); end
        if isempty(D), D = zeros(n); end
    end
    J(~isfinite(J)) = 0; D(~isfinite(D)) = 0;

    % ---------- Optional D-scaling (alpha) to tame Jeff ----------
    if ~isempty(alpha_user)
        alpha = alpha_user;
    elseif autoAlpha
        alpha = pick_alpha(J, D);
    else
        alpha = 1.0;
    end
    if alpha ~= 1, D = alpha * D; end

    % ---------- Build constant delay operator Q and Jeff ----------
    A = eye(n) - (D .* J);                 % SPM-style Hadamard delay operator base
    Acond = cond(A);
    if ~isfinite(Acond)
        Acond = Inf;
    end
    useSolve = (~isfinite(Acond) || Acond > maxCond);
    if useSolve
        Q = [];           % use A\(...) each time
        Qmode = 'solve';
    else
        % Explicit inverse is fine for n=56
        Q = inv(A); %#ok<MINV>
        Qmode = 'mul';
    end
    % Build Jeff = Q*J without forming Q if we're in solve-mode
    if strcmp(Qmode,'mul')
        Jeff = Q * J;
    elseif strcmp(Qmode,'solve')
        Jeff = A \ J;
    else
        error('Unknown Q mode: %s', Qmode);
    end

    % Spectrum diagnostics (helpful for stability)
    evJeff = eig(Jeff);
    maxRe  = max(real(evJeff));
    maxAbs = max(abs(evJeff));

    % ---------- Allocate ----------
    X = zeros(n, T);
    X(:,1) = x0;
    x      = x0;

    info.J = J; info.D = D; info.Q = Q;
    info.Qcond = condest_safe(Q, A, Qmode);
    info.Acond = Acond;
    info.method = method;
    info.notes  = 'Static-JD; supports rk4/euler/rosenbrock, substeps, and D scaling.';
    info.alpha      = alpha;
    info.evJeff     = evJeff;
    info.maxReJeff  = maxRe;
    info.maxAbsJeff = maxAbs;

    % ---------- Main loop ----------
    for k = 1:T-1
        H  = dt(k);
        tk = t(k);

        switch lower(method)
            case 'euler'
                hk = H / substeps;
                for s = 1:substeps
                    ts = tk + (s-0.5)*hk;
                    uk  = sample_u(u, ts, t, k);
                    fx  = f(x, uk, P, M);
                    fx  = apply_Q_once(fx, A, Q, Qmode);
                    if ~isempty(maxStep)
                        nx = norm(fx);
                        if nx > maxStep, fx = fx * (maxStep / nx); end
                    end
                    x = x + hk * fx;
                    if finiteCheck && any(~isfinite(x))
                        error('State became non-finite at t=%.6f (Euler)', ts);
                    end
                end

            case 'rk4'
                hk = H / substeps;
                for s = 1:substeps
                    ts1 = tk + (s-1)*hk;
                    ts2 = ts1 + 0.5*hk;
                    ts3 = ts2;
                    ts4 = ts1 + hk;

                    u1 = sample_u(u, ts1, t, k);
                    u2 = sample_u(u, ts2, t, k);
                    u3 = u2;
                    u4 = sample_u(u, ts4, t, k);

                    k1 = apply_Q_once(f(x,            u1, P, M), A, Q, Qmode);
                    k2 = apply_Q_once(f(x + 0.5*hk*k1,u2, P, M), A, Q, Qmode);
                    k3 = apply_Q_once(f(x + 0.5*hk*k2,u3, P, M), A, Q, Qmode);
                    k4 = apply_Q_once(f(x +      hk*k3,u4, P, M), A, Q, Qmode);
                    incr = (hk/6)*(k1 + 2*k2 + 2*k3 + k4);

                    if ~isempty(maxStep)
                        ni = norm(incr)/hk;
                        if ni > maxStep, incr = incr * (maxStep/ni); end
                    end
                    x = x + incr;
                    if finiteCheck && any(~isfinite(x))
                        error('State became non-finite at t=%.6f (RK4)', ts4);
                    end
                end

            case 'rosenbrock'
                % Linearly-implicit (Rosenbrock–Euler): x_{k+1} = x_k + (I - H*Jeff)^{-1} * H * (Q*f)
                hk = H / substeps;
                I  = eye(n);
                for s = 1:substeps
                    ts = tk + (s-0.5)*hk;
                    uk = sample_u(u, ts, t, k);

                    fx  = f(x, uk, P, M);
                    rhs = apply_Q_once(fx, A, Q, Qmode);   % Q*fx

                    % Solve (I - hk*Jeff)*delta = hk*rhs
                    delta = (I - hk*Jeff) \ (hk * rhs);

                    if ~isempty(maxStep)
                        nd = norm(delta)/hk;
                        if nd > maxStep, delta = delta * (maxStep/nd); end
                    end
                    x = x + delta;
                    if finiteCheck && any(~isfinite(x))
                        error('State became non-finite at t=%.6f (Rosenbrock)', ts);
                    end
                end

            otherwise
                error('Unknown method: %s', method);
        end

        X(:,k+1) = x;
    end
end

% ================= Helpers =================
function y = apply_Q_once(fx, A, Q, mode)
    switch mode
        case 'mul',   y = Q * fx;
        case 'solve', y = A \ fx;
        otherwise,    error('Unknown Q mode.');
    end
end

function uk = sample_u(u, tk, tgrid, kcol)
    if nargin < 4, kcol = []; end
    if nargin < 1 || isempty(u)
        uk = 0;
    elseif isa(u,'function_handle')
        uk = u(tk);
    elseif isstruct(u) && isfield(u,'t') && isfield(u,'U')
        tt = u.t(:).'; UU = u.U;
        if size(UU,2) ~= numel(tt), error('u.U columns must equal numel(u.t)'); end
        if isvector(UU), uk = interp1(tt, UU(:), tk, 'linear','extrap');
        else,            uk = interp1(tt, UU.', tk, 'linear','extrap').';
        end
    elseif isnumeric(u)
        if isvector(u)
            if isempty(kcol), [~,kcol] = min(abs(tgrid - tk)); end
            uk = u(min(kcol, numel(u)));
        else
            if isempty(kcol), [~,kcol] = min(abs(tgrid - tk)); end
            kcol = min(kcol, size(u,2));
            uk = u(:,kcol);
        end
    else
        error('Unsupported u format.');
    end
end

function v = getdef(s, name, def)
    if isfield(s,name), v = s.(name); else, v = def; end
end

function c = condest_safe(Q, A, mode)
    switch mode
        case 'mul',   c = cond(Q);
        case 'solve', c = cond(A);
        otherwise,    c = NaN;
    end
end

function alpha = pick_alpha(J, D)
% Back-off scaler for D to make max Re(eig(Jeff(alpha))) ≤ 0, Jeff=(I-a*D.*J)\J
    I = eye(size(J,1));
    alpha = 1.0;
    for it = 1:20
        A  = I - alpha*(D.*J);
        if any(~isfinite(A(:))) || rcond(A) < 1e-12
            alpha = alpha / 2; continue
        end
        Jeff = A \ J;
        mRe  = max(real(eig(Jeff)));
        if mRe <= 0, return; end
        alpha = alpha / 2;
    end
    % If we get here, we didn’t stabilise; return whatever we have
end
