function [J, Wsrc] = build_twocmp_observer_J_all8(M, P, weights)
% BUILD_TWOCMP_OBSERVER_J_ALL8
%   Observer for two-compartment (Vs,Vd) states in ALL 8 pops.
%
%   Returns J such that y = J * spm_vec(x)
%   - If M.L exists (nchan×ns), J = M.L * Wsrc
%   - Else J = Wsrc (one "virtual channel" per source)
%
% Model layout:
%   nk = 10: Vs=1, Vd=8 (others are conductances + STP)
%   pops: 1 SS, 2 SP, 3 SI, 4 DP, 5 DI, 6 TP, 7 RT, 8 RL/RC

    if nargin < 3 || isempty(weights)
        % --- axial "dipole-like" weights per population (default: SP+DP only)
        weights.ax = zeros(1,8);
        weights.ax(2) = 1.0;  % SP dominates
        weights.ax(4) = 0.6;  % DP secondary

        % optional small soma stabilisers (kept from your original)
        weights.bSI = 0.05; % +SI soma (small)
        weights.bDI = 0.02; % -DI soma (small)
        weights.bTP = 0.02; % +TP soma (tiny)
        weights.bRL = 0.03; % +RL soma (tiny)

        % scale axial term by gc?
        weights.use_gc = true;

        % if gc is per-pop, use it per-pop (recommended)
        weights.gc_per_pop = true;
    end

    % dims & indices
    ns = size(M.x,1); np = size(M.x,2); nk = size(M.x,3);
    iVs = 1; iVd = 8;

    % total state length after spm_vec
    N = numel(M.x);

    % helper for linear index into vec(x)
    idx = @(is, ip, ik) sub2ind([ns np nk], is, ip, ik);

    % ---- axial coupling scale(s)
    % We expect P.gc in log-space (per your priors), but allow linear too.
    if isfield(P,'gc') && weights.use_gc
        gc_raw = P.gc;
        % if it looks like log-space (common), exponentiate; otherwise keep
        % (you can force behaviour by providing P.gc already exponentiated)
        gc_lin = exp(gc_raw);
    else
        gc_lin = 1.0;
    end

    % make gc_lin either scalar or 1×8
    if isscalar(gc_lin)
        gc_pop = gc_lin * ones(1,8);
    else
        gc_lin = gc_lin(:);
        if numel(gc_lin) == 8
            gc_pop = gc_lin(:)';   % 1×8
        else
            % fallback: if weird size, just take first element
            gc_pop = gc_lin(1) * ones(1,8);
        end
    end

    % allocate triplets
    rows = []; cols = []; vals = [];

    for s = 1:ns

        % --- Axial "dipole-like" terms: sum_p ax(p) * gc(p) * (Vd - Vs)
        for p = 1:8
            a = 0;
            if isfield(weights,'ax') && numel(weights.ax) >= p
                a = weights.ax(p);
            end
            if a ~= 0
                if isfield(weights,'gc_per_pop') && weights.gc_per_pop
                    gcp = gc_pop(p);
                else
                    gcp = gc_pop(1);
                end
                rows = [rows, s, s];
                cols = [cols, idx(s,p,iVd), idx(s,p,iVs)];
                vals = [vals, gcp*a, -gcp*a];
            end
        end

        % --- Small soma stabilisers (unchanged from your original)
        if isfield(weights,'bSI') && weights.bSI ~= 0
            rows = [rows, s]; cols = [cols, idx(s,3,iVs)]; vals = [vals, +weights.bSI];
        end
        if isfield(weights,'bDI') && weights.bDI ~= 0
            rows = [rows, s]; cols = [cols, idx(s,5,iVs)]; vals = [vals, -weights.bDI];
        end
        if isfield(weights,'bTP') && weights.bTP ~= 0
            rows = [rows, s]; cols = [cols, idx(s,6,iVs)]; vals = [vals, +weights.bTP];
        end
        if isfield(weights,'bRL') && weights.bRL ~= 0
            rows = [rows, s]; cols = [cols, idx(s,8,iVs)]; vals = [vals, +weights.bRL];
        end
    end

    % ns × N sparse “per-source virtual channel” selector
    Wsrc = sparse(rows, cols, vals, ns, N);

    % Project to channels if leadfield provided
    if isfield(M,'L') && ~isempty(M.L)
        J = M.L * Wsrc;  % (nchan × N)
    else
        J = Wsrc;        % (ns × N)
    end
end
