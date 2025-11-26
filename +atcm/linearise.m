function [A,B,C,Del,x0,u0] = linearise(P, M, U)
% atcm.linearise  Linearise a thalamo–cortical model around an operating point.
%
%   [A,B,C,Del,x0,u0] = atcm.linearise(P, M, U)
%
% Expects your state-flow to be M.f(x,u,P,M) and (optionally) observation M.g(x,P,M).
% If M.f returns Jacobians, we use them:
%   [f,J,Q,D] = M.f(x,u,P,M), where J = ∂f/∂x and D may encode delays.
% Otherwise we compute A=∂f/∂x and B=∂f/∂u by central differences.
%
% Inputs
%   P : parameter struct
%   M : model struct (fields used if present)
%       .f (function handle), .g (optional), .x (operating state, any shape),
%       .u0 (baseline input), .fs (Hz), .C (obs matrix optional),
%       .obs_states (indices into vectorised x), .output_ix (row selector),
%       .Delays (struct fallback), .nu (#inputs, fallback)
%   U : optional struct with fields U.t, U.u  (baseline taken at t<=0)
%
% Outputs
%   A   : state Jacobian (nx×nx)
%   B   : input Jacobian (nx×nu)
%   C   : observation Jacobian (ny×nx)
%   Del : delays struct (if available), else empty
%   x0  : operating point (vectorised)
%   u0  : baseline input (nu×1)
%
% Notes
% - x is vectorised internally with spm_vec if available, else (:).
% - If M.g missing, C falls back to M.C, else a simple selector.
% - Designed to be stable with tiny eps and scale-aware step sizes.

% ---------- sanity ----------
if ~isfield(M,'f') || isempty(M.f)
    error('M.f (state function) is required.');
end
fhandle = M.f;

% ---------- operating point ----------
x0 = get_x0(M);
[u0, nu] = get_u0(P,M, U);

% ---------- try model-native Jacobians ----------
A = []; B = []; Del = struct();
try
    out = cell(1,4);
    [out{:}] = fhandle(x0,u0,P,M);          % [f,J,Q,D] if model supports it
    if ~isempty(out{2}), A = out{2}; end     % J = ∂f/∂x
    if ~isempty(out{4}), Del = out{4}; end   % delays struct (optional)
catch
    % model may not support multiple outputs; ignore
end

% ---------- compute Jacobians if needed ----------
if isempty(A)
    A = jacobian_x(fhandle, x0, u0, P, M);
end
if isempty(B)
    B = jacobian_u(fhandle, x0, u0, P, M, nu);
end

% ---------- observation matrix ----------
C = get_C(M, P, x0);

% ---------- fallback delays ----------
if isempty(fieldnames(Del)) && isfield(M,'Delays')
    Del = M.Delays;
end
end

% =================== helpers ===================

function x0 = get_x0(M)
% get/flatten operating point
if isfield(M,'x') && ~isempty(M.x)
    try
        x0 = spm_vec(M.x); % SPM helper if available
    catch
        x0 = M.x(:);
    end
else
    % default: 1D state with resting potential-like start
    x0 = -70;
end
x0 = x0(:);
end

function [u0, nu] = get_u0(P,M, U)
% baseline input at operating point
if nargin>=2 && ~isempty(U) && isfield(U,'t') && isfield(U,'u') && ~isempty(U.u)
    % take mean of pre-stim (t<=0); else first sample
    if any(U.t<=0)
        sel = U.t<=0;
        u0  = mean(U.u(:,sel),2);
    else
        u0  = U.u(:,1);
    end
elseif isfield(M,'u0') && ~isempty(M.u0)
    u0 = M.u0(:);
else
    % infer nu if provided
    if isfield(M,'nu') && ~isempty(M.nu)
        nu = M.nu;
    else
        nu = 1;
    end
    u0 = zeros(nu,1);
end
nu = numel(u0);
end

function A = jacobian_x(fhandle, x0, u0, P, M)
% central difference ∂f/∂x with scale-aware steps
nx = numel(x0);
fx0 = fhandle(x0, u0, P, M);
fx0 = fx0(:);
A   = zeros(numel(fx0), nx);
% step sizes
eps_abs = 1e-6;
eps_rel = 1e-4;
for i = 1:nx
    dx = zeros(nx,1);
    hi = eps_abs + eps_rel*max(1,abs(x0(i)));
    dx(i) = hi;
    fp = fhandle(x0+dx, u0, P, M);
    fm = fhandle(x0-dx, u0, P, M);
    A(:,i) = (fp(:) - fm(:)) / (2*hi);
end
end

function B = jacobian_u(fhandle, x0, u0, P, M, nu)
% central difference ∂f/∂u
fx0 = fhandle(x0, u0, P, M);
fx0 = fx0(:);
B   = zeros(numel(fx0), nu);
eps_abs = 1e-6;
eps_rel = 1e-4;
for j = 1:nu
    du = zeros(nu,1);
    hj = eps_abs + eps_rel*max(1,abs(u0(j)));
    du(j) = hj;
    fp = fhandle(x0, u0+du, P, M);
    fm = fhandle(x0, u0-du, P, M);
    B(:,j) = (fp(:) - fm(:)) / (2*hj);
end
end

function C = get_C(M, P, x0)
% build observation Jacobian C = ∂g/∂x if g exists; else use M.C or a selector
if isfield(M,'g') && ~isempty(M.g)
    try
        % try to get analytic Jacobian: [y,C] = g(x,P,M) or [y,~,C]
        out = cell(1,3);
        [out{:}] = M.g(x0, P, M);
        % try common positions
        if ~isempty(out{2}) && isnumeric(out{2})
            C = out{2};
            return
        elseif ~isempty(out{3}) && isnumeric(out{3})
            C = out{3};
            return
        end
    catch
        % numeric ∂g/∂x
        C = jacobian_g(M.g, x0, P, M);
        return
    end
end

% fallbacks
if isfield(P,'J') && ~isempty(P.J)
    C = exp(P.L)*exp(P.J);
    return
end

nx = numel(x0);
if isfield(M,'obs_states') && ~isempty(M.obs_states)
    C = zeros(numel(M.obs_states), nx);
    for k = 1:numel(M.obs_states)
        C(k, M.obs_states(k)) = 1;
    end
else
    % default: observe first state
    C = zeros(1,nx);
    C(1) = 1;
end

% optionally select a single output row
if isfield(M,'output_ix') && ~isempty(M.output_ix) && M.output_ix<=size(C,1)
    C = C(M.output_ix, :);
end
end

function C = jacobian_g(ghandle, x0, P, M)
% numeric ∂g/∂x via central difference
y0 = ghandle(x0, P, M);
y0 = y0(:);
nx = numel(x0);
C  = zeros(numel(y0), nx);
eps_abs = 1e-6;
eps_rel = 1e-4;
for i = 1:nx
    dx = zeros(nx,1);
    hi = eps_abs + eps_rel*max(1,abs(x0(i)));
    dx(i) = hi;
    yp = ghandle(x0+dx, P, M);
    ym = ghandle(x0-dx, P, M);
    C(:,i) = (yp(:) - ym(:)) / (2*hi);
end
end
