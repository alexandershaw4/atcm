function Red = tcm_auto_fix(DCMs, pE0, pC0, opts)
%TCM_AUTO_FIX  Evidence-aware greedy fixing of poorly-identified parameters.
%
% Red = tcm_auto_fix(DCMs, pE0, pC0, opts)
%
% PURPOSE
%   Given per-subject DCM fits (posterior means/covariances, free energies),
%   design a reduced model by FIXING as many parameters as possible while
%   preserving model evidence and predictive accuracy. This function:
%     1) Computes identifiability/sensitivity metrics across subjects
%     2) Ranks parameters for "fixing"
%     3) Iteratively fixes low-value parameters, guarded by evidence/error tests
%     4) Returns reduced priors (pE_red, pC_red) + a log of decisions
%
% INPUTS
%   DCMs : 1xN cell array. Each DCM{s} must include fields:
%          .Ep (posterior mean, struct OR vector matching spm_vec order)
%          .Cp (posterior covariance, full [P x P] matrix)
%          .F  (free energy / log-evidence) [optional but recommended]
%          .xY.y (observed data vector/matrix) [or supply opts.ycell]
%   pE0  : prior mean (struct or vector). If struct, spm_vec used.
%   pC0  : prior covariance (struct-with-variances or full matrix). If struct,
%          diag constructed from spm_vec(pC0).
%
%   opts : struct with optional fields
%        .simfun       : @(pStruct,DCM_s,opts) -> yhat vector (RECOMMENDED).
%                        If 'def', uses DCM.M.IS.
%        .ycell        : 1xN cell of observed data vectors if not in DCM.xY.y
%        .W            : weights for error (vector) OR 1xN cell of vectors.
%                        Defaults to ones like y.
%        .labels       : 1xP cell array param labels (optional). If absent,
%                        will generate 'theta1',...,'thetaP'.
%        .batch        : number of params to try fixing per step (default 8).
%        .tau          : allowable mean evidence drop per subject (default 1.0)
%        .eps          : allowable mean REL-MSE increase (default 0.03 = 3%).
%        .max_steps    : safety cap on greedy steps (default P).
%        .fix_value    : 'prior' | 'group_median' | 'group_mean' (default 'prior')
%        .delta_scale  : finite-difference step in SD units for sensitivity
%        .recompute    : true/false recompute scores after each accept
%        .very_small   : variance used to "fix" (you set 0)
%        .use_bmr_ig   : use IG-based analytic ΔF estimate (default true)
%        .verbose      : 0/1/2 logging (default 2 here)
%
% Extra guardrails & perf:
%        .min_free, .require_simfun, .min_IG_keep_q, .min_Sens_keep_q,
%        .critical_idx, .keep_groups, .cand_gate ('or'|'and'|'sum'),
%        .cand_frac, .min_cand
%        .use_parallel (false), .sens_topk ([] -> 30% of P), .Jfun (optional)
% Backoff:
%        .backoff_single (true), .backoff_relax_eps (0.02),
%        .backoff_relax_tau (0.2), .try_group_median_if_reject (true)
%
% OUTPUT
%   Red : struct with fields ...
%
% -------------------------------------------------------------------------

%% --------------------------- Setup & sanity ----------------------------
N = numel(DCMs);
assert(N>=1, 'DCMs must be a non-empty cell array');

has_spm = exist('spm_vec','file')==2 && exist('spm_unvec','file')==2;

if ~exist('opts','var') || isempty(opts), opts = struct; end
opts = defaulter(opts, 'batch',       2);
opts = defaulter(opts, 'tau',         1);
opts = defaulter(opts, 'eps',         0.8);
opts = defaulter(opts, 'max_steps',   Inf);
opts = defaulter(opts, 'fix_value',   'prior');%prior'
opts = defaulter(opts, 'delta_scale', 0.5);
opts = defaulter(opts, 'recompute',   true);
opts = defaulter(opts, 'very_small',  0);
opts = defaulter(opts, 'use_bmr_ig',  true);
opts = defaulter(opts, 'verbose',     2);
opts = defaulter(opts, 'simfun',     'def');

% Performance & Jacobian options
opts = defaulter(opts, 'use_parallel',  true);
opts = defaulter(opts, 'sens_topk',     []);    % if empty -> 30% of P (min 10)
opts = defaulter(opts, 'Jfun',          []);    % optional Jacobian handle

% Backoff / acceptance control
opts = defaulter(opts, 'backoff_single',             true);
opts = defaulter(opts, 'backoff_relax_eps',          0.5);
opts = defaulter(opts, 'backoff_relax_tau',          0.2);
opts = defaulter(opts, 'try_group_median_if_reject', true);
opts = defaulter(opts, 'zero_var_tol', 1e-12);   % treat <= this as fixed-by-prior

% Additional conservative guardrails (your choices kept/adapted)
opts = defaulter(opts, 'require_simfun',  true);
opts = defaulter(opts, 'min_IG_keep_q',   0.5);
opts = defaulter(opts, 'min_Sens_keep_q', 0.5);
opts = defaulter(opts, 'critical_idx',    []);
opts = defaulter(opts, 'keep_groups',     {});

opts = defaulter(opts, 'cand_gate', 'or');      % 'or' | 'and' | 'sum'
opts = defaulter(opts, 'cand_frac', 0.4);       % only for 'sum'
opts = defaulter(opts, 'min_cand',  []);        % defaults to batch
opts = defaulter(opts, 'tau_per_param', 0.5);

if isempty(pE0) && isempty(pC0)
    pE0 = DCMs{1}.M.pE;
    pC0 = DCMs{1}.M.pC;
end

% Vectorise priors
[mu0, C0, template_E, is_struct_pE, is_struct_pC] = vec_priors(pE0, pC0, has_spm);
P = numel(mu0);

opts = defaulter(opts, 'min_free',        round(sqrt(P)/2));


% Identify parameters that are already fixed by the prior (variance ~ 0)
s0_prior        = sqrt(diag(C0));               % prior SDs
fixed_by_prior  = s0_prior <= opts.zero_var_tol;
free_mask       = ~fixed_by_prior;
Kfree0          = nnz(free_mask);

% Be chatty if desired
if opts.verbose
    fprintf('[tcm_auto_fix] Ignoring %d params with prior var≈0. Free=%d/%d.\n', nnz(fixed_by_prior), Kfree0, P);
end

% Clamp min_free so we can still fix something
opts.min_free = min(opts.min_free, max(Kfree0-1, 0));


% Vectorise posteriors
[MU, Cdiag, C_full, F, ycell, Wcell, labels] = vec_posteriors(DCMs, mu0, C0, opts, has_spm);
if isfield(opts,'labels') && ~isempty(opts.labels)
    labels = opts.labels(:);
end

% Compute identifiability/sensitivity metrics
if opts.verbose, fprintf('[tcm_auto_fix] Computing parameter scores...\n'); end
%scores_tbl = compute_scores(MU, Cdiag, C_full, mu0, C0, labels, ycell, Wcell, DCMs, template_E, opts, has_spm);
scores_tbl = compute_scores(MU, Cdiag, C_full, mu0, C0, labels, ycell, Wcell, DCMs, template_E, opts, has_spm, free_mask);

% Greedy fixing
if opts.verbose
    fprintf('[tcm_auto_fix] Starting greedy fixing with tau=%.3f, eps=%.3f, batch=%d', opts.tau, opts.eps, opts.batch);
end

% If predictive checks required but simfun/y missing, abort fixing
have_sim = isfield(opts,'simfun') && ~isempty(opts.simfun) && ~isempty(ycell{1});
if opts.require_simfun && ~have_sim
    if opts.verbose
        fprintf('[tcm_auto_fix] No simfun or observed data detected. Returning FULL model (no fixing).');
    end
    Red = struct();
    Red.pE_red = unvec_param(mu0, template_E, has_spm);
    Red.pC_red = unvec_cov(C0, pC0, has_spm);
    Red.fixed_idx = []; Red.fixed_labels = {}; Red.fixed_to = [];
    Red.keep_idx = (1:P)';
    %Red.scores_initial = compute_scores(MU, Cdiag, C_full, mu0, C0, labels, ycell, Wcell, DCMs, template_E, opts, has_spm);
    Red.scores_initial = compute_scores(MU, Cdiag, C_full, mu0, C0, labels, ycell, Wcell, DCMs, template_E, opts, has_spm, free_mask);

    Red.history = table(); Red.opts = opts;
    return;
end

% Cache baseline predictions once (for predictive checks)
yhat_base = cell(1,N);
if have_sim
    for s=1:N
        yhat_base{s} = safely_sim(MU(:,s), template_E, DCMs{s}, opts, has_spm);
    end
end

keep = true(P,1);
keep(~free_mask) = false; 

fixed_vals = nan(P,1);
step = 0;
hist_step = []; hist_param_idx = {}; hist_param_lbl = {}; hist_mean_dF = []; hist_mean_dERR = []; hist_cum_fixed = []; hist_K_free = [];

cur_mu0 = mu0; cur_C0 = C0;  % running priors
MAX_STEPS = min(P, opts.max_steps);

while step < MAX_STEPS
    % Stop if we'd violate the minimum number of free parameters
    if sum(keep) <= opts.min_free, break; end

    % Next candidate pool (apply guardrails)
    cand = find(keep);
    cand = setdiff(cand, opts.critical_idx);

    ig  = scores_tbl.IG(cand);
    se  = scores_tbl.Sens(cand);
    ig_cut   = quantile(ig,   opts.min_IG_keep_q);
    sens_cut = quantile(se,   opts.min_Sens_keep_q);

    switch lower(opts.cand_gate)
        case 'and'
            mask = (ig < ig_cut) & (se < sens_cut);
        case 'sum'
            ig_n = (ig - min(ig)) ./ (max(ig)-min(ig) + eps);
            se_n = (se - min(se)) ./ (max(se)-min(se) + eps);
            combo = 0.5*ig_n + 0.5*se_n;
            thr   = quantile(combo, opts.cand_frac);
            mask  = combo <= thr;
        otherwise % 'or'
            mask = (ig < ig_cut) | (se < sens_cut);
    end

    % Relax if too few:
    min_cand = ternary(isempty(opts.min_cand), opts.batch, opts.min_cand);
    if nnz(mask) < min_cand
        %mask = (ig < ig_cut) | (se < sens_cut);
        mask = (ig <= ig_cut) | (se <= sens_cut);   % instead of < 
    end
    if nnz(mask) < min_cand
        ig_n = (ig - min(ig)) ./ (max(ig)-min(ig) + eps);
        se_n = (se - min(se)) ./ (max(se)-min(se) + eps);
        combo = 0.5*ig_n + 0.5*se_n;
        [~,ord] = sort(combo,'ascend');
        mask = false(size(cand));
        mask(ord(1:min(min_cand,numel(cand)))) = true;
    end

    cand = cand(mask);

    % Ensure each keep-group retains at least one free param
    for g = 1:numel(opts.keep_groups)
        G = intersect(opts.keep_groups{g}(:)', 1:P);
        freeG = intersect(G, find(keep));
        if numel(freeG) == 1
            cand = setdiff(cand, freeG);
        end
    end

    if isempty(cand), break; end
    [~, ord] = sort(scores_tbl.Rank(cand), 'ascend');
    pick = cand(ord(1:min(opts.batch, numel(ord))));

    % Decide fix-to values
    switch lower(opts.fix_value)
        case 'prior'
            fix_to = cur_mu0(pick);
        case 'group_median'
            fix_to = median(MU(pick,:), 2);
        case 'group_mean'
            fix_to = mean(MU(pick,:), 2);
        otherwise
            error('Unknown opts.fix_value: %s', opts.fix_value);
    end

    % Evaluate proposal (uses cached baseline, and Jfun if provided)
    [mean_dF, mean_dERR] = evaluate_fix(pick, fix_to, MU, Cdiag, mu0, cur_mu0, cur_C0, ycell, Wcell, DCMs, template_E, yhat_base, opts, has_spm);

    accept = (mean_dF >= -opts.tau) && (mean_dERR <= opts.eps);

    if opts.verbose
        fprintf(' step %d | try fixing %s | dF=%.3f dERR=%.3f -> %s\n', ...
            step+1, strjoin(labels(pick)', ', '), mean_dF, mean_dERR, ternary(accept,'ACCEPT','reject'));
    end

    % ------------------------ Backoff logic ------------------------
    if ~accept
        accepted_via_backoff = false;

        % (A) Try group-median targets if you started from 'prior'
        if opts.try_group_median_if_reject && strcmpi(opts.fix_value,'prior')
            fix_to_gm = median(MU(pick,:), 2);
            [mean_dF2, mean_dERR2] = evaluate_fix(pick, fix_to_gm, MU, Cdiag, mu0, cur_mu0, cur_C0, ycell, Wcell, DCMs, template_E, yhat_base, opts, has_spm);
            accept2 = (mean_dF2 >= -(opts.tau + opts.backoff_relax_tau)) && (mean_dERR2 <= (opts.eps + opts.backoff_relax_eps));
            if opts.verbose
                fprintf('   backoff A: group-median target | dF=%.3f dERR=%.3f -> %s\n', mean_dF2, mean_dERR2, ternary(accept2,'ACCEPT','reject'));
            end
            if accept2
                fix_to = fix_to_gm; mean_dF = mean_dF2; mean_dERR = mean_dERR2; accept = true; accepted_via_backoff = true;
            end
        end

        % (B) Try best single parameter from the batch
        if ~accept && opts.backoff_single
            best = struct('idx',[], 'fix',[], 'dF',-Inf, 'dERR',Inf, 'ok',false);
            for j = 1:numel(pick)
                pj   = pick(j);
                fixj = fix_to(j);
                [dFj, dERRj] = evaluate_fix(pj, fixj, MU, Cdiag, mu0, cur_mu0, cur_C0, ycell, Wcell, DCMs, template_E, yhat_base, opts, has_spm);
                okj = (dFj >= -(opts.tau + opts.backoff_relax_tau)) && (dERRj <= (opts.eps + opts.backoff_relax_eps));
                better = (dERRj < best.dERR - 1e-12) || ((abs(dERRj - best.dERR) <= 1e-12) && (dFj > best.dF));
                if better
                    best = struct('idx', pj, 'fix', fixj, 'dF', dFj, 'dERR', dERRj, 'ok', okj);
                end
            end
            if opts.verbose && ~isempty(best.idx)
                fprintf('   backoff B: best single %s | dF=%.3f dERR=%.3f -> %s\n', labels{best.idx}, best.dF, best.dERR, ternary(best.ok,'ACCEPT','reject'));
            end
            if best.ok
                pick   = best.idx;
                fix_to = best.fix;
                mean_dF = best.dF; mean_dERR = best.dERR; accept = true; accepted_via_backoff = true;
            end
        end

        if ~accepted_via_backoff
            break; % stopping rule remains
        end
    end

    % Accept: update priors and bookkeeping
    step = step + 1;
    [cur_mu0, cur_C0] = apply_fix_to_priors(cur_mu0, cur_C0, pick, fix_to, opts.very_small);
    keep(pick) = false; fixed_vals(pick) = fix_to;

    % Log
    hist_step(end+1,1) = step; %#ok<AGROW>
    hist_param_idx{end+1,1} = pick; %#ok<AGROW>
    hist_param_lbl{end+1,1} = labels(pick); %#ok<AGROW>
    hist_mean_dF(end+1,1)   = mean_dF; %#ok<AGROW>
    hist_mean_dERR(end+1,1) = mean_dERR; %#ok<AGROW>
    hist_cum_fixed(end+1,1) = sum(~keep); %#ok<AGROW>
    hist_K_free(end+1,1)    = sum(keep); %#ok<AGROW>

    % Optionally recompute scores under current priors (posteriors unchanged)
    if opts.recompute
        scores_tbl = compute_scores(MU, Cdiag, C_full, cur_mu0, cur_C0, labels, ycell, Wcell, DCMs, template_E, opts, has_spm, free_mask);
    end

    % Stop if we'd drop below the minimum number of free parameters next step
    if sum(keep) <= opts.min_free, break; end
end

% Prepare outputs
Red = struct();
Red.pE_red = unvec_param(cur_mu0, template_E, has_spm);
Red.pC_red = unvec_cov(cur_C0, pC0, has_spm);
Red.fixed_idx    = find(~keep);
Red.fixed_labels = labels(~keep);
Red.fixed_to     = fixed_vals(~keep);
Red.keep_idx     = find(keep);
Red.scores_initial = scores_tbl;

if isempty(hist_step)
    Red.history = table();
else
    Red.history = table(hist_step, hist_param_idx, hist_param_lbl, hist_mean_dF, hist_mean_dERR, hist_cum_fixed, hist_K_free, ...
        'VariableNames', {'step','param_idx','param_labels','mean_dF','mean_dERR','cum_fixed','K_free'});
end

Red.opts = opts;

if opts.verbose && ~isempty(hist_step)
    fprintf('[tcm_auto_fix] Finished. Fixed %d/%d params.\n', numel(Red.fixed_idx), P);
end

end % main function

%% ============================ SUBFUNCTIONS ==============================
function opts = defaulter(opts, field, val)
    if ~isfield(opts, field) || isempty(opts.(field))
        opts.(field) = val;
    end
end

function t = ternary(cond, a, b)
    if cond, t = a; else, t = b; end
end

function [mu0, C0, template_E, is_struct_pE, is_struct_pC] = vec_priors(pE0, pC0, has_spm)
    is_struct_pE = isstruct(pE0);
    is_struct_pC = isstruct(pC0);
    if has_spm
        if is_struct_pE
            mu0 = spm_vec(pE0);
            template_E = pE0;
        else
            mu0 = pE0(:);
            template_E = [];
        end
        if is_struct_pC
            v = spm_vec(pC0); C0 = diag(v(:));
        else
            C0 = pC0;
        end
    else
        % Fallback: flatten by field order (must be consistent everywhere!)
        if is_struct_pE
            [mu0, template_E] = fallback_vec_struct(pE0);
        else
            mu0 = pE0(:); template_E = [];
        end
        if is_struct_pC
            v = fallback_vec_struct(pC0); C0 = diag(v(:));
        else
            C0 = pC0;
        end
    end
end

function [MU, Cdiag, C_full, F, ycell, Wcell, labels] = vec_posteriors(DCMs, mu0, C0, opts, has_spm)
    N = numel(DCMs); P = numel(mu0);
    MU = zeros(P,N); Cdiag = zeros(P,N); C_full = cell(N,1); F = nan(N,1);

    ycell = cell(1,N); Wcell = cell(1,N);

    % labels fall-back
    labels = arrayfun(@(i) sprintf('theta%d', i), 1:P, 'uni', false)';

    for s = 1:N
        Ep = DCMs{s}.Ep; Cp = DCMs{s}.Cp;

        % if the stored Cp is only for active parameters, re-embed in full space
        if length(Cp) ~= length(Cdiag)
            Cp = atcm.fun.reembedreducedcovariancematrix(DCMs{s},Cp);
        end

        if has_spm && isstruct(Ep)
            mu = spm_vec(Ep);
        else
            mu = Ep(:);
        end
        MU(:,s) = mu(:);
        C_full{s} = Cp;
        Cdiag(:,s) = diag(Cp);
        if isfield(DCMs{s},'F'), F(s) = DCMs{s}.F; end

        % Observations
        if isfield(DCMs{s},'xY') && isfield(DCMs{s}.xY,'y') && ~isempty(DCMs{s}.xY.y)
            ycell{s} = DCMs{s}.xY.y;
        elseif isfield(opts,'ycell') && ~isempty(opts.ycell)
            ycell{s} = opts.ycell{s};
        else
            ycell{s} = [];
        end
        if isfield(opts,'W') && ~isempty(opts.W)
            if iscell(opts.W), Wcell{s} = opts.W{s}; else, Wcell{s} = opts.W; end
        else
            if isempty(ycell{s}), Wcell{s} = []; else, Wcell{s} = ones(numel(ycell{s}),1); end
        end
    end

    % attempt to derive labels from pE structure if available
    if has_spm && isstruct(DCMs{1}.Ep)
        try
            [~, labels] = spm_vec(DCMs{1}.Ep);
        catch
        end
    end
end

function scores_tbl = compute_scores(MU, Cdiag, C_full, mu0, C0, labels, ycell, Wcell, DCMs, template_E, opts, has_spm, free_mask)
    [P,N] = size(MU);
    if nargin < 13 || isempty(free_mask), free_mask = true(P,1); end
    free_idx  = find(free_mask);
    fixed_idx = find(~free_mask);

    % Preallocate (full-length); compute on free only
    IG       = nan(P,1);
    Delta    = nan(P,1);
    Rho      = nan(P,1);
    Sens_med = nan(P,1);
    BwOverW  = nan(P,1);
    Rank     = inf(P,1);                 % fixed-by-prior get Rank = +Inf so they’re ignored

    % --- Core scores on free subset ---
    s0 = sqrt(diag(C0) + eps);
    IG_free    = 0.5 * mean(log((s0(free_idx).^2) ./ (Cdiag(free_idx,:) + eps)), 2);
    Delta_free = mean(abs(MU(free_idx,:) - mu0(free_idx)) ./ (s0(free_idx) + eps), 2);

    % Posterior collinearity (compute full, then pick free rows)
    Rho_full = zeros(P,1);
    for s=1:N
        R = corrcov(C_full{s} + 1e-12*eye(P));
        R(1:P+1:end) = 0;
        Rho_full = Rho_full + max(abs(R),[],2) / N;
    end
    Rho_free = Rho_full(free_idx);

    % Sensitivity
    Sens = zeros(P,N);
    have_sim = isfield(opts,'simfun') && ~isempty(opts.simfun);
    if ~isempty(opts.Jfun) && ~isempty(ycell{1})
        for s=1:N
            y = ycell{s}; if isempty(y), continue; end
            w = Wcell{s}; if isempty(w), w = ones(numel(y),1); end
            mu = MU(:,s);
            pStruct = unvec_param(mu, template_E, has_spm);
            J = opts.Jfun(pStruct, DCMs{s}, opts);     % [numel(y) x P]
            WJ = J .* w(:);
            Sens(:,s) = sqrt(sum(WJ.^2,1))';
        end
    elseif have_sim && ~isempty(ycell{1})
        % Focus FD on weakest free params
        IGn = normalize(IG_free,   'range');
        Dn  = normalize(Delta_free,'range');
        baseRank_free = (1-IGn) + (1-Dn);
        if isempty(opts.sens_topk)
            sens_topk = max(10, round(1*numel(free_idx)));
        else
            sens_topk = opts.sens_topk;
        end
        [~,ordFD] = sort(baseRank_free, 'descend');
        sens_idx = free_idx(ordFD(1:min(sens_topk, numel(free_idx))));

        parfor s = 1:N
            Sens_col = zeros(P,1);

            y = ycell{s};
            if isempty(y)
                Sens(:,s) = Sens_col; continue;
            end
            w = Wcell{s}; if isempty(w), w = ones(numel(y),1); end

            mu  = MU(:,s);
            y1  = safely_sim(mu, template_E, DCMs{s}, opts, has_spm);

            for ii = 1:numel(sens_idx)
                p    = sens_idx(ii);
                mu2  = mu;
                step = opts.delta_scale * s0(p);
                if step == 0, continue; end   % guard even if mask missed something
                mu2(p) = mu2(p) + step;

                y2 = safely_sim(mu2, template_E, DCMs{s}, opts, has_spm);
                Sens_col(p) = norm((y2(:) - y1(:)) .* w(:));
            end

            Sens(:,s) = Sens_col;
        end
    end
    Sens_med_free = median(Sens(free_idx,:), 2);

    % Between/within variance (free only)
    Bw_free      = var(MU(free_idx,:), 0, 2);
    Ww_free      = median(Cdiag(free_idx,:), 2);
    BwOverW_free = Bw_free ./ (Ww_free + 1e-12);

    % Normalize each score (free rows only)
    S_free = [IG_free, Delta_free, Rho_free, Sens_med_free, BwOverW_free];
    S_free = normalize(S_free, 1, 'range');

    Rank_free = (1 - S_free(:,1)) + (1 - S_free(:,2)) + S_free(:,3) + (1 - S_free(:,4)) + (1 - S_free(:,5));

    % Scatter back into full-length vectors
    IG(free_idx)       = IG_free;
    Delta(free_idx)    = Delta_free;
    Rho(free_idx)      = Rho_free;
    Sens_med(free_idx) = Sens_med_free;
    BwOverW(free_idx)  = BwOverW_free;
    Rank(free_idx)     = Rank_free;

    % Build table (flag fixed-by-prior)
    fixed_by_prior = ~free_mask;
    scores_tbl = table((1:P)', labels, IG, Delta, Rho, Sens_med, BwOverW, Rank, fixed_by_prior, ...
        'VariableNames', {'idx','param','IG','Delta','Rho','Sens','BwOverW','Rank','fixed_by_prior'});
    scores_tbl = sortrows(scores_tbl, 'Rank', 'ascend');
end


function yhat = safely_sim(mu_vec, template_E, DCMs_s, opts, has_spm)
    try
        pStruct = unvec_param(mu_vec, template_E, has_spm);
        if ischar(opts.simfun) && strcmp(opts.simfun,'def')
            yhat = spm_vec(DCMs_s.M.IS(pStruct, DCMs_s.M, DCMs_s.xU));
        else
            yhat = opts.simfun(pStruct, DCMs_s, opts);
        end
        if isempty(yhat), yhat = zeros(1,1); end
    catch ME
        warning('simfun failed: %s. Returning zeros.', ME.message);
        yhat = zeros(1,1);
    end
end

function [mean_dF, mean_dERR] = evaluate_fix(pick, fix_to, MU, Cdiag, mu0, cur_mu0, cur_C0, ycell, Wcell, DCMs, template_E, yhat_base, opts, has_spm)
    [P,N] = size(MU);

    % Evidence change: IG-based proxy (sum over picked)
    dF = zeros(N,1);
    if opts.use_bmr_ig
        s0 = sqrt(diag(cur_C0)+eps);
        IG_all = 0.5 * log((s0.^2) ./ (Cdiag + eps)); % P x N
        dF = - mean(sum(IG_all(pick,:),1)', 1);
    else
        dF(:) = 0;
    end

    % % Evidence change: IG-based proxy (average per fixed param)
    % dF = zeros(N,1);
    % if opts.use_bmr_ig
    %     s0 = sqrt(diag(cur_C0) + eps);
    %     IG_all  = 0.5 * log((s0.^2) ./ (Cdiag + eps));    % P x N
    %     IG_pick = IG_all(pick,:);                          % |pick| x N
    %     dF = -mean(IG_pick, 1)';                           % average over picked params, per subject
    % else
    %     dF(:) = 0;
    % end

    % Predictive error change: no refit; either linearized (Jfun) or one sim per subject
    have_sim = isfield(opts,'simfun') && ~isempty(opts.simfun) && ~isempty(ycell{1});
    dERR = zeros(N,1);
    % ----- Predictive error change (robust, relative to baseline with floor) -----
    have_sim = isfield(opts,'simfun') && ~isempty(opts.simfun) && ~isempty(ycell{1});
    dERR = zeros(N,1);
    if have_sim
        for s=1:N
            y = ycell{s}; if isempty(y), continue; end
            w = Wcell{s}; if isempty(w), w = ones(numel(y),1); end
            mu = MU(:,s);

            yhat_old = yhat_base{s};  % cached baseline

            if ~isempty(opts.Jfun)
                pStruct = unvec_param(mu, template_E, has_spm);
                J  = opts.Jfun(pStruct, DCMs{s}, opts);         % [numel(y) x P]
                dt = zeros(numel(mu),1); dt(pick) = fix_to - mu(pick);
                yhat_new = yhat_old(:) + J*dt;
            else
                mu_new   = mu; mu_new(pick) = fix_to;
                yhat_new = safely_sim(mu_new, template_E, DCMs{s}, opts, has_spm);
            end

            % --- length guard ---
            yo = y(:); yh_old = yhat_old(:); yh_new = yhat_new(:);
            M = min([numel(yo), numel(yh_old), numel(yh_new), numel(w)]);
            yo = yo(1:M); yh_old = yh_old(1:M); yh_new = yh_new(1:M); ww = w(1:M); ww = ww(:);

            yh_old = spm_vec(yh_old);
            yo = spm_vec(yo);
            
            % --- weighted MSEs ---
            e_old = mean(((yh_old - yo).^2) .* ww);
            e_new = mean(((yh_new - yo).^2) .* ww);

            % --- robust denominator: floor at tiny fraction of signal power ---
            sig_pow = mean((yo.^2) .* ww);                      % scale proxy
            e_floor = max(1e-6 * sig_pow, 1e-12);               % floor
            denom   = max(e_old, e_floor);

            dERR(s) = (e_new - e_old) / denom;                  % fractional increase
            % clip pathological blow-ups to keep a single subject from vetoing everything
            dERR(s) = max(min(dERR(s), 5), -0.9);
        end
    else
        dERR(:) = 0;
    end

    mean_dF = mean(dF);
    mean_dERR = mean(dERR);
end

function e = rel_mse(yhat, y, w, yhat_ref)
    % relative MSE; if yhat_ref given, normalise by its error
    yhat = spm_vec(yhat(:));
    y    = spm_vec(y(:));

    if numel(yhat) ~= numel(y)
        M = min(numel(yhat), numel(y));
        yhat = yhat(1:M); y = y(1:M); w = w(1:M);
    end
    if nargin<3 || isempty(w), w = ones(numel(y),1); end
    w = w(:);

    num = mean(((yhat - y).^2) .* w);

    if nargin >= 4 && ~isempty(yhat_ref)
        den = mean(((yhat_ref - y).^2) .* w) + 1e-12;   % baseline error
    else
        den = mean((y.^2) .* w) + 1e-12;                 % fallback: signal power
    end
    e = num / den;
end

function [mu0_new, C0_new] = apply_fix_to_priors(mu0, C0, pick, fix_to, very_small)
    mu0_new = mu0; mu0_new(pick) = fix_to(:);
    C0_new  = C0;
    C0_new(pick,:) = 0; C0_new(:,pick) = 0; 
    for i = 1:numel(pick)
        C0_new(pick(i), pick(i)) = very_small;
    end
end

function S = unvec_param(mu_vec, template_E, has_spm)
    if has_spm && ~isempty(template_E)
        S = spm_unvec(mu_vec, template_E);
    else
        S = mu_vec; % best-effort fallback
    end
end

function C = unvec_cov(C0_vec, pC0_template, has_spm)
    if isstruct(pC0_template)
        if has_spm
            v = diag(C0_vec);
            C = spm_unvec(v, pC0_template); % returns struct of variances
        else
            C = C0_vec; % fallback: keep numeric
        end
    else
        C = C0_vec; % numeric stays numeric
    end
end

function [v, template] = fallback_vec_struct(S)
    % VERY simplified fallback: concatenates numeric fields in order
    f = fieldnames(S);
    v = [];
    for i=1:numel(f)
        x = S.(f{i});
        v = [v; x(:)]; %#ok<AGROW>
    end
    template = S;
end

%% =============================== USAGE =================================
% Example (pseudo-code):
%
% opts = struct();
% opts.simfun = @(pStruct,DCM_s,opts) Alex_LaplaceTFwD(pStruct, DCM_s); % your forward
% opts.batch = 2; opts.tau = 1.0; opts.eps = 0.03; opts.fix_value = 'prior';
% Red = tcm_auto_fix(DCMs, DCMs{1}.M.pE, DCMs{1}.M.pC, opts);
%
% Outputs:
%   Red.pE_red, Red.pC_red  -> use as DCM.M.pE / DCM.M.pC for re-fitting
%   Red.history             -> per-step log (dF proxy and dERR)
%   Red.scores_initial      -> parameter ranks & diagnostics
%
% Notes:
% - Evidence proxy uses the drop in information gain when fixing params.
%   If you have a proper BMR routine, you can replace evaluate_fix() with a
%   call that computes analytic ΔF using (μ,Σ) and prior changes per Friston
%   et al., 2016 (Bayesian Model Reduction). The plug-in is designed so you
%   can swap it with a warm-start VL re-fit as well.
% - Sensitivity uses finite differences with step = delta_scale * prior SD.
% - If SPM isn’t on the path, vectorisation falls back to a simplistic
%   concatenation; ensure consistent ordering across priors and posteriors.
%
% -------------------------------------------------------------------------
