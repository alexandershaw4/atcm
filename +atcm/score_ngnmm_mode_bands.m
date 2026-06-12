function [score, info] = score_ngnmm_mode_bands(P, M, targetBands, opts)
% score_ngnmm_mode_bands
% -------------------------------------------------------------------------
% Score whether an NGNMM/DCM operating regime has useful, stable, resonant,
% and distinct modes across requested frequency bands.
%
% Intended use:
%   [score, info] = score_ngnmm_mode_bands(P, M, targetBands, opts)
%
% Lower score is better.
%
% Inputs
% -------------------------------------------------------------------------
% P           parameter structure, e.g. DCM.M.pE or candidate P
% M           model structure with fields M.x and M.f
% targetBands nbands x 2 matrix, e.g. [2 5; 6 10; 10 14; 15 25; 30 45]
% opts        optional structure
%
% Useful opts fields
% -------------------------------------------------------------------------
% opts.fpTol              fixed-point residual tolerance, default 1e-4
% opts.fpSearchTol         tolerance passed to find_fixed_point_robust, default 1e-8
% opts.minQ               desired minimum Q for matched modes, default 0.75
% opts.maxQ               upper soft limit for Q, default 8
% opts.minFreq            lower allowed mode frequency, default 0.5 Hz
% opts.maxFreq            upper allowed mode frequency, default 90 Hz
% opts.requireDistinct     use each eigenmode once only, default true
% opts.highFreqPenalty     add explicit penalty for missing high bands, default true
% opts.verbose            print diagnostics, default false
%
% Outputs
% -------------------------------------------------------------------------
% score       scalar objective value; lower is better
% info        diagnostic structure with eigenvalues, frequencies, Q, matches
%
% AS / ChatGPT drop-in, 2026

    if nargin < 3 || isempty(targetBands)
        targetBands = [ ...
             2   5
             6  10
            10  14
            15  25
            30  45
            45  65];
    end

    if nargin < 4 || isempty(opts)
        opts = struct();
    end

    opts = local_defaults(opts);

    info = struct();
    info.failed = true;
    info.score = Inf;
    info.targetBands = targetBands;

    % ---------------------------------------------------------------------
    % Fixed point / operating point
    % ---------------------------------------------------------------------
    if ~isfield(M,'x') || isempty(M.x)
        score = 1e12;
        info.reason = 'M.x missing';
        return;
    end

    x_template = M.x;

    try
        x_fp = atcm.fun.find_fixed_point_robust(P, M, opts.fpSearchTol);
        x_fp = spm_unvec(x_fp(:), x_template);
    catch ME
        score = 1e12;
        info.reason = ['fixed point failed: ' ME.message];
        return;
    end

    try
        [f,A,D] = feval(M.f, x_fp, 0, P, M); %#ok<ASGLU>
    catch ME
        score = 1e12;
        info.reason = ['model evaluation failed: ' ME.message];
        return;
    end

    fp_resid = norm(f(:));
    info.fp_resid = fp_resid;
    info.x_fp = x_fp;

    if ~isfinite(fp_resid) || fp_resid > opts.fpTol
        score = 1e9 + 1e6*safe_num(fp_resid);
        info.reason = 'bad fixed point residual';
        info.score = score;
        return;
    end

    % ---------------------------------------------------------------------
    % Eigenmode diagnostics
    % ---------------------------------------------------------------------
    ev = eig(A);
    freq = abs(imag(ev)) ./ (2*pi);
    damp = real(ev);
    Q = abs(imag(ev)) ./ max(eps, -2*damp);

    info.ev = ev;
    info.freq = freq;
    info.damp = damp;
    info.Q = Q;

    if any(~isfinite(ev)) || any(~isfinite(freq)) || any(~isfinite(Q))
        score = 1e12;
        info.reason = 'non-finite eigen diagnostics';
        info.score = score;
        return;
    end

    % Stability penalty: all real parts should be negative.
    unstable = real(ev) > -opts.stabilityMargin;
    if any(unstable)
        stability_penalty = opts.stabilityWeight * sum((real(ev(unstable)) + opts.stabilityMargin).^2 + 1);
    else
        stability_penalty = 0;
    end

    % Positive-frequency representatives only, to avoid selecting both sides
    % of a conjugate pair.
    ix = find(imag(ev) > 0 & ...
              freq >= opts.minFreq & freq <= opts.maxFreq & ...
              real(ev) < -opts.stabilityMargin & ...
              Q > opts.absoluteMinQ);

    if isempty(ix)
        score = 1e8 + stability_penalty + 1e3*fp_resid;
        info.reason = 'no usable positive-frequency modes';
        info.score = score;
        return;
    end

    f_use = freq(ix);
    q_use = Q(ix);
    r_use = real(ev(ix)); %#ok<NASGU>

    nb = size(targetBands,1);
    matchedLocal = nan(nb,1);
    matchedIdx = nan(nb,1);
    matchedHz = nan(nb,1);
    matchedQ = nan(nb,1);
    matchedReal = nan(nb,1);
    matchedInsideBand = false(nb,1);

    used = false(numel(ix),1);

    band_cost = 0;
    missing_penalty = 0;
    distinct_penalty = 0;
    q_cost = 0;

    % ---------------------------------------------------------------------
    % Match one distinct mode to each band
    % ---------------------------------------------------------------------
    for b = 1:nb
        lo = targetBands(b,1);
        hi = targetBands(b,2);
        centre = sqrt(lo*hi);

        if opts.requireDistinct
            candidates = find(~used);
        else
            candidates = (1:numel(ix))';
        end

        if isempty(candidates)
            missing_penalty = missing_penalty + opts.missingBandPenalty;
            continue;
        end

        in_band = candidates(f_use(candidates) >= lo & f_use(candidates) <= hi);

        if ~isempty(in_band)
            cand = in_band;
            dist = abs(log(f_use(cand) + eps) - log(centre));
            qpen = max(0, opts.minQ - q_use(cand)).^2;
            overq = max(0, q_use(cand) - opts.maxQ).^2;

            local_cost = opts.freqWeight .* dist.^2 + ...
                         opts.qUnderWeight .* qpen + ...
                         opts.qOverWeight .* overq;

            [bestCost,ii] = min(local_cost);
            chosen = cand(ii);
            band_cost = band_cost + bestCost;
            matchedInsideBand(b) = true;
        else
            cand = candidates;
            dist_to_band = zeros(size(cand));

            for cc = 1:numel(cand)
                ff = f_use(cand(cc));
                if ff < lo
                    dist_to_band(cc) = abs(log(ff + eps) - log(lo));
                elseif ff > hi
                    dist_to_band(cc) = abs(log(ff + eps) - log(hi));
                else
                    dist_to_band(cc) = 0;
                end
            end

            qpen = max(0, opts.minQ - q_use(cand)).^2;
            overq = max(0, q_use(cand) - opts.maxQ).^2;

            local_cost = opts.outOfBandWeight .* dist_to_band.^2 + ...
                         opts.qUnderWeight .* qpen + ...
                         opts.qOverWeight .* overq + ...
                         opts.outOfBandPenalty;

            [bestCost,ii] = min(local_cost);
            chosen = cand(ii);
            band_cost = band_cost + bestCost;
            missing_penalty = missing_penalty + opts.outOfBandPenalty;
            matchedInsideBand(b) = false;
        end

        if opts.requireDistinct
            used(chosen) = true;
        end

        matchedLocal(b) = chosen;
        matchedIdx(b) = ix(chosen);
        matchedHz(b) = f_use(chosen);
        matchedQ(b) = q_use(chosen);
        matchedReal(b) = real(ev(ix(chosen)));
    end

    % Penalise duplicate matches if distinctness has been disabled.
    if ~opts.requireDistinct
        valid = matchedIdx(~isnan(matchedIdx));
        distinct_penalty = opts.duplicatePenalty * (numel(valid) - numel(unique(valid)));
    end

    % Explicit Q cost for matched modes.
    for b = 1:nb
        if ~isnan(matchedQ(b))
            q = matchedQ(b);
            if q < opts.minQ
                q_cost = q_cost + opts.qUnderWeight * (opts.minQ - q)^2;
            elseif q > opts.maxQ
                q_cost = q_cost + opts.qOverWeight * (q - opts.maxQ)^2;
            end
        else
            q_cost = q_cost + opts.missingBandPenalty;
        end
    end

    % ---------------------------------------------------------------------
    % Optional high-frequency coverage penalty
    % ---------------------------------------------------------------------
    hf_penalty = 0;
    if opts.highFreqPenalty
        has_30_45 = any(f_use >= 30 & f_use <= 45 & q_use >= opts.highFreqMinQ);
        has_45_65 = any(f_use >= 45 & f_use <= 65 & q_use >= opts.gammaMinQ);

        if any(targetBands(:,1) <= 45 & targetBands(:,2) >= 30) && ~has_30_45
            hf_penalty = hf_penalty + opts.highFreqMissingPenalty;
        end

        if any(targetBands(:,1) <= 65 & targetBands(:,2) >= 45) && ~has_45_65
            hf_penalty = hf_penalty + opts.highFreqMissingPenalty;
        end
    end

    % Softly discourage all matched modes collapsing into narrow range.
    spread_penalty = 0;
    validHz = matchedHz(~isnan(matchedHz));
    if numel(validHz) >= 2
        logSpread = range(log(validHz + eps));
        desiredSpread = range(log(sqrt(prod(targetBands,2)) + eps));
        spread_penalty = opts.spreadWeight * max(0, 0.5*desiredSpread - logSpread)^2;
    end

    score = band_cost + q_cost + missing_penalty + distinct_penalty + ...
            hf_penalty + spread_penalty + stability_penalty + ...
            opts.fpWeight * fp_resid;

    % ---------------------------------------------------------------------
    % Package diagnostics
    % ---------------------------------------------------------------------
    info.failed = false;
    info.reason = 'ok';
    info.score = score;
    info.stability_penalty = stability_penalty;
    info.band_cost = band_cost;
    info.q_cost = q_cost;
    info.missing_penalty = missing_penalty;
    info.hf_penalty = hf_penalty;
    info.spread_penalty = spread_penalty;
    info.usableIdx = ix;
    info.usableHz = f_use;
    info.usableQ = q_use;
    info.matchedLocal = matchedLocal;
    info.matchedIdx = matchedIdx;
    info.matchedHz = matchedHz;
    info.matchedQ = matchedQ;
    info.matchedReal = matchedReal;
    info.matchedInsideBand = matchedInsideBand;

    if opts.verbose
        fprintf('score %.4f | fp %.3e\n', score, fp_resid);
        Tmatch = table(targetBands(:,1), targetBands(:,2), matchedHz, matchedQ, matchedInsideBand, ...
            'VariableNames', {'BandLow','BandHigh','MatchedHz','MatchedQ','InsideBand'});
        disp(Tmatch);
    end
end

% =========================================================================
% Defaults
% =========================================================================
function opts = local_defaults(opts)
    if ~isfield(opts,'fpTol'), opts.fpTol = 1e-4; end
    if ~isfield(opts,'fpSearchTol'), opts.fpSearchTol = 1e-8; end
    if ~isfield(opts,'fpWeight'), opts.fpWeight = 1e3; end

    if ~isfield(opts,'minQ'), opts.minQ = 0.75; end
    if ~isfield(opts,'absoluteMinQ'), opts.absoluteMinQ = 0.03; end
    if ~isfield(opts,'maxQ'), opts.maxQ = 8; end
    if ~isfield(opts,'highFreqMinQ'), opts.highFreqMinQ = 0.50; end
    if ~isfield(opts,'gammaMinQ'), opts.gammaMinQ = 0.30; end

    if ~isfield(opts,'minFreq'), opts.minFreq = 0.5; end
    if ~isfield(opts,'maxFreq'), opts.maxFreq = 90; end

    if ~isfield(opts,'freqWeight'), opts.freqWeight = 1; end
    if ~isfield(opts,'outOfBandWeight'), opts.outOfBandWeight = 6; end
    if ~isfield(opts,'outOfBandPenalty'), opts.outOfBandPenalty = 2; end
    if ~isfield(opts,'missingBandPenalty'), opts.missingBandPenalty = 10; end
    if ~isfield(opts,'duplicatePenalty'), opts.duplicatePenalty = 5; end

    if ~isfield(opts,'qUnderWeight'), opts.qUnderWeight = 3; end
    if ~isfield(opts,'qOverWeight'), opts.qOverWeight = 0.2; end

    if ~isfield(opts,'highFreqPenalty'), opts.highFreqPenalty = true; end
    if ~isfield(opts,'highFreqMissingPenalty'), opts.highFreqMissingPenalty = 8; end

    if ~isfield(opts,'spreadWeight'), opts.spreadWeight = 2; end

    if ~isfield(opts,'stabilityMargin'), opts.stabilityMargin = 1e-6; end
    if ~isfield(opts,'stabilityWeight'), opts.stabilityWeight = 1e6; end

    if ~isfield(opts,'requireDistinct'), opts.requireDistinct = true; end
    if ~isfield(opts,'verbose'), opts.verbose = false; end
end

function x = safe_num(x)
    if isempty(x) || ~isfinite(x)
        x = 1e6;
    end
end
