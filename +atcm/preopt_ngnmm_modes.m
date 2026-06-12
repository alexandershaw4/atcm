function [Pbest, infoBest, hist] = preopt_ngnmm_modes(P0, M, targetBands, opts)
% preopt_ngnmm_modes
% -------------------------------------------------------------------------
% Pre-optimise an NGNMM prior parameter set so that the fixed-point
% linearisation has stable, resonant, and distinct modes across specified
% frequency bands.
%
% Intended to run before aFitDCM/VL.
%
% Example
% -------------------------------------------------------------------------
% targetBands = [2 5; 6 10; 10 14; 15 25; 30 45; 45 65];
% opts = struct('maxIter',120,'nStarts',8,'display',true);
% [pE_mode, modeInfo, hist] = preopt_ngnmm_modes(DCM.M.pE, DCM.M, targetBands, opts);
% DCM.M.pE = pE_mode;
%
% Notes
% -------------------------------------------------------------------------
% This does not fit data amplitude.  It only tunes the prior dynamical
% regime so VL begins from a model with useful modes.
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

    if ~isfield(M,'Npop')
        error('preopt_ngnmm_modes:MissingNpop', 'M.Npop is required.');
    end

    Np = M.Npop;

    if isfield(M,'Kmask')
        Kmask = M.Kmask;
    else
        Kmask = ones(Np,Np);
    end

    kidx = find(Kmask(:) == 1);

    % Ensure required parameter fields exist.
    P0 = local_complete_P(P0, M, Kmask);

    % Pack parameters selected for pre-optimisation.
    spec = local_make_spec(P0, M, Kmask, opts);
    theta0 = local_pack(P0, spec);

    % Clamp initial vector.
    theta0 = local_clamp_theta(theta0, spec, opts);

    % Objective function.
    obj = @(theta) local_objective(theta, P0, M, targetBands, opts, spec);

    hist = struct();
    hist.theta0 = theta0;
    hist.targetBands = targetBands;
    hist.starts = [];
    hist.scores = [];
    hist.exitflag = [];
    hist.output = {};

    bestScore = Inf;
    thetaBest = theta0;
    infoBest = struct('failed',true,'reason','not evaluated');

    % Build starting points. First is deterministic P0, rest jittered.
    starts = zeros(numel(theta0), opts.nStarts);
    starts(:,1) = theta0;

    rngState = rng;
    if ~isempty(opts.rngSeed)
        rng(opts.rngSeed);
    end

    for s = 2:opts.nStarts
        starts(:,s) = local_jitter_theta(theta0, spec, opts);
    end

    if ~isempty(opts.rngSeed)
        rng(rngState);
    end

    % ---------------------------------------------------------------------
    % Optimise from each start
    % ---------------------------------------------------------------------
    if opts.display
        fprintf('\n--------------- NGNMM MODE PRE-OPTIMISATION ---------------\n');
        fprintf('Starts: %d | maxIter/start: %d\n', opts.nStarts, opts.maxIter);
        disp(array2table(targetBands, 'VariableNames', {'BandLow','BandHigh'}));
    end

    for s = 1:opts.nStarts
        th0 = starts(:,s);

        if opts.display
            fprintf('\n[start %d/%d]\n', s, opts.nStarts);
        end

        switch lower(opts.method)
            case 'fminsearch'
                options = optimset( ...
                    'Display','off', ...
                    'MaxIter',opts.maxIter, ...
                    'MaxFunEvals',opts.maxFunEvals, ...
                    'TolX',opts.tolX, ...
                    'TolFun',opts.tolFun);

                [th, sc, exitflag, output] = fminsearch(obj, th0, options);

            case 'pattern'
                % Very simple coordinate/pattern search fallback without toolbox.
                [th, sc, exitflag, output] = local_pattern_search(obj, th0, opts);

            otherwise
                error('Unknown optimisation method: %s', opts.method);
        end

        th = local_clamp_theta(th, spec, opts);
        P = local_unpack(th, P0, spec, opts);
        [sc2, info] = atcm.score_ngnmm_mode_bands(P, M, targetBands, opts.scoreOpts);

        hist.starts(s).theta0 = th0; %#ok<AGROW>
        hist.starts(s).theta = th; %#ok<AGROW>
        hist.starts(s).score = sc2; %#ok<AGROW>
        hist.starts(s).rawScore = sc; %#ok<AGROW>
        hist.starts(s).info = info; %#ok<AGROW>
        hist.exitflag(s) = exitflag; %#ok<AGROW>
        hist.output{s} = output; %#ok<AGROW>
        hist.scores(s) = sc2; %#ok<AGROW>

        if opts.display
            fprintf('score %.6f | fp %.3e | reason: %s\n', ...
                sc2, getfield_default(info,'fp_resid',NaN), getfield_default(info,'reason',''));
            local_print_match(info);
        end

        if sc2 < bestScore
            bestScore = sc2;
            thetaBest = th;
            infoBest = info;
        end
    end

    Pbest = local_unpack(thetaBest, P0, spec, opts);
    [~, infoBest] = atcm.score_ngnmm_mode_bands(Pbest, M, targetBands, opts.scoreOpts);
    infoBest.thetaBest = thetaBest;
    infoBest.spec = spec;
    infoBest.bestScore = bestScore;

    if opts.display
        fprintf('\nMode pre-optimisation complete\n');
        fprintf('Best score: %.6f\n', bestScore);
        fprintf('Fixed-point residual: %.3e\n', infoBest.fp_resid);
        local_print_match(infoBest);
        local_print_modes(infoBest);
    end
end

% =========================================================================
% Parameter packing specification
% =========================================================================
function spec = local_make_spec(P0, M, Kmask, opts)
    Np = M.Npop;

    spec = struct();
    spec.Np = Np;
    spec.include = struct();
    spec.names = {};
    spec.slices = struct();

    n = 0;

    if opts.fitLogQIFRate
        spec.slices.logQIFRate = n + (1:1);
        spec.names{end+1} = 'logQIFRate';
        n = n + 1;
    end

    if opts.fitLogAlpha
        spec.slices.logAlpha = n + (1:Np);
        spec.names{end+1} = 'logAlpha';
        n = n + Np;
    end

    if opts.fitEta
        spec.slices.eta = n + (1:Np);
        spec.names{end+1} = 'eta';
        n = n + Np;
    end

    if opts.fitLogDelta
        spec.slices.logDelta = n + (1:Np);
        spec.names{end+1} = 'logDelta';
        n = n + Np;
    end

    if opts.fitLogKappa
        kidx = find(Kmask(:) == 1);
        spec.kidx = kidx;
        spec.slices.logKappa = n + (1:numel(kidx));
        spec.names{end+1} = 'logKappa';
        n = n + numel(kidx);
    else
        spec.kidx = [];
    end

    spec.nTheta = n;

    % Keep original dimensions/templates.
    spec.Ptemplate = P0;
end

function theta = local_pack(P, spec)
    theta = zeros(spec.nTheta,1);

    if isfield(spec.slices,'logQIFRate')
        theta(spec.slices.logQIFRate) = P.logQIFRate;
    end

    if isfield(spec.slices,'logAlpha')
        theta(spec.slices.logAlpha) = P.logAlpha(:);
    end

    if isfield(spec.slices,'eta')
        theta(spec.slices.eta) = P.eta(:);
    end

    if isfield(spec.slices,'logDelta')
        theta(spec.slices.logDelta) = P.logDelta(:);
    end

    if isfield(spec.slices,'logKappa')
        LK = P.logKappa;
        theta(spec.slices.logKappa) = LK(spec.kidx);
    end
end

function P = local_unpack(theta, P0, spec, opts)
    P = P0;
    theta = local_clamp_theta(theta, spec, opts);

    if isfield(spec.slices,'logQIFRate')
        P.logQIFRate = theta(spec.slices.logQIFRate);
    end

    if isfield(spec.slices,'logAlpha')
        P.logAlpha = reshape(theta(spec.slices.logAlpha), 1, spec.Np);
    end

    if isfield(spec.slices,'eta')
        P.eta = reshape(theta(spec.slices.eta), 1, spec.Np);
    end

    if isfield(spec.slices,'logDelta')
        P.logDelta = reshape(theta(spec.slices.logDelta), 1, spec.Np);
    end

    if isfield(spec.slices,'logKappa')
        LK = P.logKappa;
        LK(spec.kidx) = theta(spec.slices.logKappa);
        P.logKappa = LK;
    end
end

function theta = local_clamp_theta(theta, spec, opts)
    theta = theta(:);

    if isfield(spec.slices,'logQIFRate')
        sl = spec.slices.logQIFRate;
        theta(sl) = min(max(theta(sl), log(opts.bounds.qifRate(1))), log(opts.bounds.qifRate(2)));
    end

    if isfield(spec.slices,'logAlpha')
        sl = spec.slices.logAlpha;
        theta(sl) = min(max(theta(sl), log(opts.bounds.alpha(1))), log(opts.bounds.alpha(2)));
    end

    if isfield(spec.slices,'eta')
        sl = spec.slices.eta;
        theta(sl) = min(max(theta(sl), opts.bounds.eta(1)), opts.bounds.eta(2));
    end

    if isfield(spec.slices,'logDelta')
        sl = spec.slices.logDelta;
        theta(sl) = min(max(theta(sl), log(opts.bounds.delta(1))), log(opts.bounds.delta(2)));
    end

    if isfield(spec.slices,'logKappa')
        sl = spec.slices.logKappa;
        theta(sl) = min(max(theta(sl), log(opts.bounds.kappa(1))), log(opts.bounds.kappa(2)));
    end
end

function theta = local_jitter_theta(theta0, spec, opts)
    theta = theta0;

    if isfield(spec.slices,'logQIFRate')
        sl = spec.slices.logQIFRate;
        theta(sl) = theta(sl) + opts.jitter.logQIFRate * randn(numel(sl),1);
    end

    if isfield(spec.slices,'logAlpha')
        sl = spec.slices.logAlpha;
        theta(sl) = theta(sl) + opts.jitter.logAlpha * randn(numel(sl),1);
    end

    if isfield(spec.slices,'eta')
        sl = spec.slices.eta;
        theta(sl) = theta(sl) + opts.jitter.eta * randn(numel(sl),1);
    end

    if isfield(spec.slices,'logDelta')
        sl = spec.slices.logDelta;
        theta(sl) = theta(sl) + opts.jitter.logDelta * randn(numel(sl),1);
    end

    if isfield(spec.slices,'logKappa')
        sl = spec.slices.logKappa;
        theta(sl) = theta(sl) + opts.jitter.logKappa * randn(numel(sl),1);
    end

    theta = local_clamp_theta(theta, spec, opts);
end

function val = local_objective(theta, P0, M, targetBands, opts, spec)
    try
        P = local_unpack(theta, P0, spec, opts);
        [val, info] = atcm.score_ngnmm_mode_bands(P, M, targetBands, opts.scoreOpts); %#ok<ASGLU>
    catch
        val = 1e12;
    end

    if ~isfinite(val)
        val = 1e12;
    end
end

% =========================================================================
% Simple pattern search fallback
% =========================================================================
function [theta, best, exitflag, output] = local_pattern_search(obj, theta0, opts)
    theta = theta0(:);
    best = obj(theta);

    step = opts.patternInitialStep;
    nEval = 1;
    it = 0;

    while it < opts.maxIter && step > opts.patternMinStep && nEval < opts.maxFunEvals
        it = it + 1;
        improved = false;

        for j = 1:numel(theta)
            for dir = [-1 1]
                th = theta;
                th(j) = th(j) + dir * step;
                sc = obj(th);
                nEval = nEval + 1;

                if sc < best
                    best = sc;
                    theta = th;
                    improved = true;
                end

                if nEval >= opts.maxFunEvals
                    break;
                end
            end
            if nEval >= opts.maxFunEvals
                break;
            end
        end

        if ~improved
            step = step * 0.5;
        end
    end

    exitflag = double(step <= opts.patternMinStep);
    output = struct('iterations',it,'funcCount',nEval,'finalStep',step);
end

% =========================================================================
% Defaults / helpers
% =========================================================================
function opts = local_defaults(opts)
    if ~isfield(opts,'method'), opts.method = 'fminsearch'; end
    if ~isfield(opts,'maxIter'), opts.maxIter = 120; end
    if ~isfield(opts,'maxFunEvals'), opts.maxFunEvals = 10000; end
    if ~isfield(opts,'tolX'), opts.tolX = 1e-3; end
    if ~isfield(opts,'tolFun'), opts.tolFun = 1e-3; end
    if ~isfield(opts,'nStarts'), opts.nStarts = 6; end
    if ~isfield(opts,'display'), opts.display = true; end
    if ~isfield(opts,'rngSeed'), opts.rngSeed = []; end

    if ~isfield(opts,'fitLogQIFRate'), opts.fitLogQIFRate = true; end
    if ~isfield(opts,'fitLogAlpha'), opts.fitLogAlpha = true; end
    if ~isfield(opts,'fitEta'), opts.fitEta = true; end
    if ~isfield(opts,'fitLogDelta'), opts.fitLogDelta = true; end
    if ~isfield(opts,'fitLogKappa'), opts.fitLogKappa = true; end

    if ~isfield(opts,'bounds') || isempty(opts.bounds)
        opts.bounds = struct();
    end
    if ~isfield(opts.bounds,'qifRate'), opts.bounds.qifRate = [30 250]; end
    if ~isfield(opts.bounds,'alpha'), opts.bounds.alpha = [10 140]; end
    if ~isfield(opts.bounds,'eta'), opts.bounds.eta = [-4 3]; end
    if ~isfield(opts.bounds,'delta'), opts.bounds.delta = [0.02 1.0]; end
    if ~isfield(opts.bounds,'kappa'), opts.bounds.kappa = [1e-6 0.8]; end

    if ~isfield(opts,'jitter') || isempty(opts.jitter)
        opts.jitter = struct();
    end
    if ~isfield(opts.jitter,'logQIFRate'), opts.jitter.logQIFRate = 0.35; end
    if ~isfield(opts.jitter,'logAlpha'), opts.jitter.logAlpha = 0.35; end
    if ~isfield(opts.jitter,'eta'), opts.jitter.eta = 0.5; end
    if ~isfield(opts.jitter,'logDelta'), opts.jitter.logDelta = 0.25; end
    if ~isfield(opts.jitter,'logKappa'), opts.jitter.logKappa = 0.45; end

    if ~isfield(opts,'patternInitialStep'), opts.patternInitialStep = 0.25; end
    if ~isfield(opts,'patternMinStep'), opts.patternMinStep = 1e-3; end

    % Score options are passed through to score_ngnmm_mode_bands.
    if ~isfield(opts,'scoreOpts') || isempty(opts.scoreOpts)
        opts.scoreOpts = struct();
    end

    % Copy relevant top-level display setting into score options if absent.
    if ~isfield(opts.scoreOpts,'verbose'), opts.scoreOpts.verbose = false; end
end

function P = local_complete_P(P, M, Kmask)
    Np = M.Npop;

    if ~isfield(P,'logQIFRate'), P.logQIFRate = log(75); end
    if ~isfield(P,'logAlpha'), P.logAlpha = log(40) * ones(1,Np); end
    if ~isfield(P,'eta'), P.eta = zeros(1,Np); end
    if ~isfield(P,'logDelta'), P.logDelta = log(0.2) * ones(1,Np); end

    if ~isfield(P,'logKappa')
        P.logKappa = log(1e-6 * ones(Np,Np));
        P.logKappa(Kmask == 1) = log(0.05);
    end
end

function local_print_match(info)
    if ~isfield(info,'targetBands') || ~isfield(info,'matchedHz')
        return;
    end

    Tmatch = table( ...
        info.targetBands(:,1), ...
        info.targetBands(:,2), ...
        info.matchedHz(:), ...
        info.matchedQ(:), ...
        info.matchedInsideBand(:), ...
        'VariableNames', {'BandLow','BandHigh','MatchedHz','MatchedQ','InsideBand'});
    disp(Tmatch);
end

function local_print_modes(info)
    if ~isfield(info,'ev')
        return;
    end

    ev = info.ev;
    freq = abs(imag(ev)) ./ (2*pi);
    damp = real(ev);
    Q = abs(imag(ev)) ./ max(eps, -2*damp);

    T = table(real(ev), imag(ev), freq, Q, ...
        'VariableNames', {'Real','Imag','Hz','Q'});

    ix = find(imag(ev) > 0 & freq > 0.5 & freq < 90 & real(ev) < 0 & Q > 0.03);
    if ~isempty(ix)
        fprintf('\nUsable positive-frequency modes, sorted by Q:\n');
        disp(sortrows(T(ix,:), 'Q', 'descend'));
    end
end

function out = getfield_default(s, field, default)
    if isstruct(s) && isfield(s,field)
        out = s.(field);
    else
        out = default;
    end
end
