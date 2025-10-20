function [y, x, t, info] = dcm_int_delay(P, M, U)
% DCM-compatible time-domain integrator with state-state delays
%
% [y, x, t, info] = dcm_int_delay(P, M, U)
%
% PURPOSE
%   Numerically integrates DCM state equations in the time domain while
%   handling pairwise state delays returned as the 3rd output of M.f.
%   This enables fitting transient phenomena (e.g., sleep spindles) where
%   conduction/synaptic delays matter.
%
% INTERFACES (SPM/DCM style)
%   M.f    - handle to state-evolution function
%            [fx, dfdx, D] = M.f(x, u, P, M)
%              x    : n×1 current state vector
%              u    : m×1 current input
%              P    : parameter structure or vector
%              M    : model struct
%              fx   : n×1 state derivative (dx/dt)
%              dfdx : (optional) Jacobian; ignored by this integrator
%              D    : n×n non-negative delay matrix (seconds)
%                     D(i,j) = delay from state j to i
%                     May be empty, scalar (uniform), or full matrix.
%
%   M.g    - (optional) observation function
%            y = M.g(x, P, M) -> p×1 observation
%
%   M.x    - (optional) n×1 initial state; default zeros(n,1)
%   M.dt   - (optional) integrator step in seconds; default 1/1000 s
%   M.t    - (optional) vector of sample times (sec) for outputs; if
%            omitted, constructed from U.dt and size(U.u,2)
%   M.nonlinear_delay_interp - (optional) 'linear' (default) | 'cubic'
%   M.integrator             - (optional) 'rk4' (default) | 'euler' | 'heun'
%   M.state_clip             - (optional) [min max] to clip states each step
%   M.progress               - (optional) true/false; print progress
%
%   U.u    - m×T matrix of inputs sampled at U.dt or with U.t provided
%   U.dt   - scalar seconds between input samples (if U.t absent)
%   U.t    - 1×T time vector (sec) for inputs (optional)
%
% RETURNS
%   y    : p×Ns samples (observations). If M.g absent, y = x (identity)
%   x    : n×Ns state trajectory at sample times t
%   t    : 1×Ns vector of output/sample times (sec)
%   info : struct with metadata (dt_int, method, max_delay, n_steps,...)
%
% NOTES
% • Delays are handled via a method-of-steps with a rolling history buffer
%   and (by default) linear interpolation of past states. The delay matrix
%   may vary over time if returned by M.f; it is re-evaluated each step.
% • The integrator sub-steps at M.dt and records at the times in t.
% • To seed a non-zero past for delays, set M.x_hist (n×K) and M.t_hist
%   (1×K) covering at least max(D) seconds before the first time in t.
%
% EXAMPLE
%   % minimal usage with constructed output grid from U
%   [y, x, t] = dcm_int_delay(P, M, U);
%
% -------------------------------------------------------------------------

% ---- Basic checks
assert(isfield(M,'f') && isa(M.f,'function_handle'), 'M.f must be a function handle');

% ---- Output time grid
if isfield(M,'t') && ~isempty(M.t)
    tout = M.t(:)';
elseif isfield(U,'t') && ~isempty(U.t)
    tout = U.t(:)';
elseif isfield(U,'dt') && ~isempty(U.dt) && isfield(U,'u') && ~isempty(U.u)
    Nt = size(U.u,2);
    tout = (0:Nt-1)*U.dt;
else
    %error('Provide M.t or U.t or U.dt with U.u to define output/sample times.');
end

if isfield(U,'dt') && isfield(M,'pst')
     tout = M.pst(:)';
     mu   = 1 * exp(P.d(2));
     sig  = 1.5+exp(P.d(3));
     U.u  = exp(-0.5*((tout-mu)/sig).^2);
end

% ---- Integrator configuration
if ~isfield(M,'dt') || isempty(M.dt)
    M.dt = 1/1000; % 1 ms default substep
end
method = getfield_with_default(M,'integrator','rk4');
interp = getfield_with_default(M,'nonlinear_delay_interp','linear');
clip   = getfield_with_default(M,'state_clip',[]);
progress = isfield(M,'progress') && M.progress;

% ---- Initial state
% First call to M.f to infer dimension if needed
n = infer_state_dim(M, P);
if ~isfield(M,'x') || isempty(M.x)
    x0 = zeros(n,1);
else
    x0 = M.x(:);
    assert(numel(x0)==n, 'M.x dimension mismatch with state size');
end

% ---- Inputs as function of time
u_fun = make_input_function(U);

% ---- Pre-build history buffer
[t_hist, X_hist] = seed_history(tout(1), x0, M, P, n);

% ---- Step through time using method-of-steps
Nt_out = numel(tout);
X_out  = zeros(n, Nt_out);
Y_out  = [];

% if observation function exists, preflight to get p
has_g = isfield(M,'g') && isa(M.g,'function_handle');
if has_g
    try p = numel(M.g(x0, P, M)); catch, p = n; end
    Y_out = zeros(p, Nt_out);
end

h = M.dt;                 % integrator substep
current_t = t_hist(end);  % may start < tout(1) if seeded
k_out = 1;

% main integration loop until we reach the last requested output time
n_steps = 0; max_delay_seen = 0;

while k_out <= Nt_out
    target_t = tout(k_out);
    % integrate in sub-steps until we reach target_t (or are at it already)
    while current_t + 1e-12 < target_t
        % Current input
        u = u_fun(current_t);

        % Build delayed state vector x_d given D returned by M.f at x
        % First, get instantaneous D to know how far back we need
        [~, ~, D] = safe_call_f(M.f, X_hist(:,end), u, P, M);
        D = normalize_delay_matrix(D, n);
        maxD = max(D(:));
        max_delay_seen = max(max_delay_seen, maxD);

        % Construct delayed composite state where each i may draw from a
        % different time for each presynaptic j. We achieve this by
        % building x_tau for the RHS evaluation via the helper below.
        get_x_tau = @(x_now) delayed_state(x_now, D, current_t, t_hist, X_hist, interp);

        % One sub-step with chosen method
        switch lower(method)
            case 'euler'
                xt = X_hist(:,end);
                fx = rhs_f(M.f, xt, u, P, M, get_x_tau);
                x_next = xt + h*fx;
            case 'heun' % improved Euler
                xt = X_hist(:,end);
                k1 = rhs_f(M.f, xt, u, P, M, get_x_tau);
                x_pred = xt + h*k1;
                u_mid = u_fun(current_t + h);
                k2 = rhs_f(M.f, x_pred, u_mid, P, M, get_x_tau);
                x_next = xt + 0.5*h*(k1 + k2);
            otherwise % 'rk4'
                xt = X_hist(:,end);
                k1 = rhs_f(M.f, xt, u, P, M, get_x_tau);
                u2 = u_fun(current_t + 0.5*h);
                k2 = rhs_f(M.f, xt + 0.5*h*k1, u2, P, M, get_x_tau);
                k3 = rhs_f(M.f, xt + 0.5*h*k2, u2, P, M, get_x_tau);
                u4 = u_fun(current_t + h);
                k4 = rhs_f(M.f, xt + h*k3, u4, P, M, get_x_tau);
                x_next = xt + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
        end

        % optional clipping for stability
        if ~isempty(clip)
            x_next = min(max(x_next, clip(1)), clip(2));
        end

        % advance time and append to history buffer
        current_t = current_t + h;
        t_hist = [t_hist, current_t]; %#ok<AGROW>
        X_hist = [X_hist, x_next];     %#ok<AGROW>
        n_steps = n_steps + 1;
    end

    % We have arrived at (or surpassed by eps) the recording time
    % Use linear interpolation over the last two history samples
    X_out(:,k_out) = interpolate_state(t_hist, X_hist, tout(k_out), interp);

    if has_g
        Y_out(:,k_out) = M.g(X_out(:,k_out), P, M);
    end

    if progress && (mod(k_out, max(1, floor(Nt_out/20)))==0)
        fprintf('dcm\\_int\\_delay: %3.0f%% (t=%.3f s)\n', 100*k_out/Nt_out, tout(k_out));
    end

    k_out = k_out + 1;
end

% ---- Pack outputs
x = X_out;
if has_g, y = Y_out; else, y = X_out; end

t = tout;
info = struct('dt_int',h,'method',method,'interp',interp, ...
              'max_delay_seen',max_delay_seen,'n_steps',n_steps);

end

% ======================================================================
% Helpers
% ======================================================================
function n = infer_state_dim(M, P)
    % Try to infer n by calling M.f at zeros with zero input
    if isfield(M,'n') && ~isempty(M.n)
        n = M.n; return; end
    try
        xz = zeros(1,1); % try 1 first
        u0 = zeros(input_dim(M),1);
        fx = M.f(xz, u0, P, M);
        n = numel(fx);
    catch
        % fallback: try n from M.x
        if isfield(M,'x') && ~isempty(M.x)
            n = numel(M.x);
        else
            error('Cannot infer state dimension; please set M.n or M.x.');
        end
    end
end

function m = input_dim(M)
    if isfield(M,'m') && ~isempty(M.m)
        m = M.m; return; end
    m = 1; % default single input
end

function val = getfield_with_default(S, name, default)
    if isfield(S, name) && ~isempty(S.(name))
        val = S.(name);
    else
        val = default;
    end
end

function [t_hist, X_hist] = seed_history(t0, x0, M, P, n)
    % Prepare initial history; use provided M.x_hist/t_hist if present,
    % otherwise create a single point at t0 with x0.
    if isfield(M,'x_hist') && ~isempty(M.x_hist) && isfield(M,'t_hist') && ~isempty(M.t_hist)
        t_hist = M.t_hist(:)';
        X_hist = M.x_hist;
        assert(size(X_hist,1)==n, 'M.x_hist must be n×K');
        assert(numel(t_hist)==size(X_hist,2), 'M.t_hist and M.x_hist length mismatch');
        % ensure last time == t0; if not, extend with steady integration
        if abs(t_hist(end) - t0) > 10*eps
            if t_hist(end) > t0
                error('Provided history ends after first output time');
            end
            % pad with final state up to t0
            t_hist = [t_hist, t0]; %#ok<AGROW>
            X_hist = [X_hist, X_hist(:,end)]; %#ok<AGROW>
        end
    else
        t_hist = t0;
        X_hist = x0(:);
    end
end

function x_tau = delayed_state(x_now, D, t_now, t_hist, X_hist, interp)
    % Build the delayed state used by RHS evaluation.
    % For additive coupling models of the form dx/dt = f(x_tau, ...), a
    % common approximation is to evaluate all presynaptic contributions at
    % their respective delayed times. Here we assemble an effective state
    % vector x_tau where entry j is the state j evaluated at t_now - d_j,
    % using the maximum incoming delay to j (max over i of D(i,j)). This
    % preserves each channel's own delay footprint while keeping vector
    % dimension n. If you require i-by-j specific sampling, embed that
    % inside M.f using access to the full D and an interpolation helper.
    n = size(D,1);
    x_tau = zeros(n,1);
    % use column-wise (per presynaptic j) max incoming delay
    dj = max(D,[],1)'; % n×1
    for j = 1:n
        tj = t_now - max(dj(j),0);
        x_tau(j) = interpolate_state(t_hist, X_hist, tj, interp, j);
    end
end

function v = interpolate_state(t_hist, X_hist, t_query, interp, j_idx)
    % Interpolate full state vector (or one index) at t_query.
    if nargin < 5 || isempty(j_idx)
        % vector output
        v = zeros(size(X_hist,1),1);
        for k = 1:size(X_hist,1)
            v(k) = interp1_local(t_hist, X_hist(k,:), t_query, interp);
        end
    else
        v = interp1_local(t_hist, X_hist(j_idx,:), t_query, interp);
    end
end

function val = interp1_local(t_hist, x_hist_row, tq, method)
    if tq <= t_hist(1)
        val = x_hist_row(1);
        return
    end
    if tq >= t_hist(end)
        val = x_hist_row(end);
        return
    end
    switch lower(method)
        case 'cubic'
            val = interp1(t_hist, x_hist_row, tq, 'pchip');
        otherwise % 'linear'
            val = interp1(t_hist, x_hist_row, tq, 'linear');
    end
end

function fx = rhs_f(fhandle, x_now, u, P, M, get_x_tau)
    % Evaluate RHS using delayed state proxy
    x_tau = get_x_tau(x_now); %#ok<NASGU> 
    % Convention: pass the delayed vector via M.x_tau if the model wishes
    % to use it; otherwise, many models define f using x directly.
    if ~isfield(M,'use_x_tau') || M.use_x_tau
        M_local = M; %#ok<NASGU>
        try
            fx = fhandle(x_tau, u, P, M); % prefer delayed x
        catch
            % fall back to current state if model signature disallows
            fx = fhandle(x_now, u, P, M);
        end
    else
        fx = fhandle(x_now, u, P, M);
    end
    fx = fx(:);
end

function D = normalize_delay_matrix(Din, n)
    if isempty(Din)
        D = zeros(n);
    elseif isscalar(Din)
        D = repmat(max(Din,0), n, n);
    else
        D = max(Din,0);
        if ~isequal(size(D), [n n])
            error('Delay matrix D must be n×n');
        end
    end
end

function [fx, dfdx, D] = safe_call_f(fhandle, x, u, P, M)
    try
        [fx, dfdx, D] = fhandle(x, u, P, M);
    catch
        % allow fewer outputs
        outs = cell(1,3);
        [outs{:}] = fhandle(x, u, P, M);
        fx   = outs{1};
        dfdx = []; if numel(outs)>=2, dfdx = outs{2}; end
        D    = []; if numel(outs)>=3, D    = outs{3}; end
    end
    if isempty(D), D = zeros(numel(x)); end
end

function u_fun = make_input_function(U)
    if nargin<1 || isempty(U) || ~isfield(U,'u') || isempty(U.u)
        u_fun = @(t) 0; % scalar zero input
        return
    end
    Uu = U.u; % m×T
    if isfield(U,'t') && ~isempty(U.t)
        Ut = U.t(:)';
    elseif isfield(U,'dt') && ~isempty(U.dt)
        Ut = (0:size(Uu,2)-1)*U.dt;
    else
        error('Provide U.t or U.dt with U.u');
    end
    % vectorised interpolation over time for each input channel
    u_fun = @(tq) interp_inputs(Ut, Uu, tq);
end

function u = interp_inputs(Ut, Uu, tq)
    m = size(Uu,1);
    u = zeros(m,1);
    for i = 1:m
        if tq <= Ut(1)
            u(i) = Uu(i,1);
        elseif tq >= Ut(end)
            u(i) = Uu(i,end);
        else
            u(i) = interp1(Ut, Uu(i,:), tq, 'linear');
        end
    end
end
