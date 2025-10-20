function [t, X, Y, Meta] = integrate_tcm_modular(P, M, dt, T, ufun, x0, obs)
%INTEGRATE_TCM_MODULAR  Fixed-step RK4 integrator for the extended TCM.
%
% [t, X, Y, Meta] = integrate_tcm_modular(P, M, dt, T, ufun, x0, obs)
%
% Inputs
%   P     : parameter struct (fields matching M.labels, log-params allowed)
%   M     : model struct with fields:
%             .f    -> function handle: [f,J,D] = tcm_modular(xvec, u, P, M)
%             .Reg  -> registries (Reg.Pops, Reg.Channels, etc.)
%             .x    -> (optional) state template sized [ns, np, 1+numel(channels)]
%   dt    : time step in seconds (e.g. 1e-3)
%   T     : total duration in seconds (e.g. 10)
%   ufun  : function handle u = ufun(t)  (column vector); use @(t)[0;0] if none
%   x0    : (optional) initial state vectorised with spm_vec; if empty, zeros
%   obs   : (optional) observation struct with fields:
%             .C -> observation row vector on spm_vec(x) (defaults to mean of voltages)
%
% Outputs
%   t     : time vector (1 x N)
%   X     : state trajectory [N x Nx] (each row is spm_vec of state at time t(n))
%   Y     : observation vector [N x 1] (if obs provided/auto)
%   Meta  : struct with useful info (sizes, index maps)
%
% Notes
% - Assumes delays (if any) are handled internally by tcm_modular using M.Reg/M.Comp.
% - RK4 is used for a good stability/speed compromise. You can swap to Euler/Heun if needed.
%
% Dr Alexander D. Shaw — extended TCM integrator

% --- Sizes / state template ---
%if ~isfield(M,'Reg') || ~isfield(M,'f')
%    error('M must contain .Reg and .f (@tcm_modular).');
%end
%Reg   = M.Reg;
[ns,np,nk] = size(M.x);

%ns    = 1;                        % your current code uses one source (ns=1)
%np    = numel(Reg.Pops);
%nk    = 1 + numel(Reg.Channels);  % [V, gates...]
Nx    = ns * np * nk;             % total state dimension (vectorised)

% State template
if ~isfield(M,'x') || isempty(M.x)
    M.x = zeros(ns, np, nk);
end

% Initial condition
if nargin < 6 || isempty(x0)
    x = spm_vec(M.x);             % zeros by default
else
    x = x0(:);
    if numel(x) ~= Nx
        error('x0 has %d elements but expected %d.', numel(x), Nx);
    end
end

% Time & storage
N  = max(1, round(T/dt));
t  = (0:N-1) * dt;
X  = zeros(N, Nx);
Y  = [];
compute_y = (nargin >= 7) && ~isempty(obs);
if compute_y
    if ~isfield(obs,'C') || isempty(obs.C)
        % Default observation: mean of voltages across populations (state #1 in each pop)
        C = zeros(1, Nx);
        for p = 1:np
            idx = sub2ind([ns, np, nk], 1, p, 1); %#ok<NASGU>  % (not used directly)
            % spm_vec order is linear; use helper below to get index of V for pop p
        end
        % Build a mask picking the first state of each population (V)
        Vmask = false(1, Nx);
        for p = 1:np
            Vmask( vec_index(ns, np, nk, 1, p, 1) ) = true;
        end
        C = double(Vmask) / sum(Vmask); % mean of voltages
        obs.C = C;
    end
    Y = zeros(N, 1);
end

% --- RK4 integrator ---
for n = 1:N
    tn = t(n);
    u  = ufun(tn);
    k1 = f_eval(x, u, P, M);

    u2 = ufun(tn + 0.5*dt);
    k2 = f_eval(x + 0.5*dt*k1, u2, P, M);

    k3 = f_eval(x + 0.5*dt*k2, u2, P, M);

    u4 = ufun(tn + dt);
    k4 = f_eval(x + dt*k3, u4, P, M);

    x  = x + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
    X(n,:) = x(:).';

    if compute_y
        Y(n,1) = obs.C * x;
    end
end

% Meta / indices helper
Meta = struct();
Meta.ns = ns; Meta.np = np; Meta.nk = nk; Meta.Nx = Nx;
Meta.index = @(s,p,k) vec_index(ns, np, nk, s, p, k);

end

% ===== helper: single f evaluation =====
function dx = f_eval(xvec, u, P, M)
% tcm_modular returns [f,J,D], we need just f
[f, ~, ~] = M.f(xvec, u, P, M);
dx = f;
end

% ===== helper: index in spm_vec order for (source s, population p, state k) =====
function ii = vec_index(ns, np, nk, s, p, k)
% Linear index assuming memory order matches spm_vec on [ns x np x nk]
% i.e., kk runs fastest, then p, then s.
ii = (s-1)*(np*nk) + (p-1)*nk + k;
end
