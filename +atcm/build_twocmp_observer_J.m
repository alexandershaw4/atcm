function [J, Wsrc] = build_twocmp_observer_J(M, P, weights)
% BUILD_TWOCMP_OBSERVER_J
%   Returns J such that y = J * spm_vec(x)
%   - If M.L is present (channels×sources), J = M.L * Wsrc
%   - Else J = Wsrc (one "virtual" channel per source)
%
% Layout matches your model:
%   states per pop (nk=10): Vs=1, Vd=8
%   populations: 1 SS, 2 SP, 3 SI, 4 DP, 5 DI, 6 TP, 7 RT, 8 RL

    if nargin < 3 || isempty(weights)
        weights.aSP = 1.0;  % SP dipole weight
        weights.aDP = 0.6;  % DP dipole weight
        weights.bSI = 0.05; % +SI soma (small)
        weights.bDI = 0.02; % -DI soma (small)
        weights.bTP = 0.02; % +TP soma (tiny)
        weights.bRL = 0.03; % +RL soma (tiny)
        weights.use_gc = true; % include axial coupling scale in J
    end

    % dims & indices
    ns = size(M.x,1); np = size(M.x,2); nk = size(M.x,3);
    iVs = 1; iVd = 8;

    % axial coupling scale
    if isfield(P,'gc') && weights.use_gc
        gc = exp(P.gc);
    else
        gc = 1.0;
    end

    % total state length after spm_vec
    N = numel(M.x);

    % row i = source i (one virtual channel per source)
    rows = []; cols = []; vals = [];

    % helper for linear index into vec(x)
    idx = @(is, ip, ik) sub2ind([ns np nk], is, ip, ik);

    for s = 1:ns
        % SP (pop 2): gc*(Vd - Vs)
        rows = [rows, s, s];
        cols = [cols, idx(s,2,iVd), idx(s,2,iVs)];
        vals = [vals, gc*weights.aSP, -gc*weights.aSP];

        % DP (pop 4): gc*(Vd - Vs)
        rows = [rows, s, s];
        cols = [cols, idx(s,4,iVd), idx(s,4,iVs)];
        vals = [vals, gc*weights.aDP, -gc*weights.aDP];

        % Small soma stabilisers (signs as suggested)
        if weights.bSI ~= 0
            rows = [rows, s]; cols = [cols, idx(s,3,iVs)]; vals = [vals,  +weights.bSI];
        end
        if weights.bDI ~= 0
            rows = [rows, s]; cols = [cols, idx(s,5,iVs)]; vals = [vals,  -weights.bDI];
        end
        if weights.bTP ~= 0
            rows = [rows, s]; cols = [cols, idx(s,6,iVs)]; vals = [vals,  +weights.bTP];
        end
        if weights.bRL ~= 0
            rows = [rows, s]; cols = [cols, idx(s,8,iVs)]; vals = [vals,  +weights.bRL];
        end
    end

    % ns × N sparse “per-source virtual channel” selector
    Wsrc = sparse(rows, cols, vals, ns, N);

    % Project to channels if leadfield provided
    if isfield(M,'L') && ~isempty(M.L)
        % M.L: (nchan × ns)
        J = M.L * Wsrc;  % -> (nchan × N)
    else
        J = Wsrc;        % -> (ns × N) one virtual channel per source
    end
end
