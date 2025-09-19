function Pri = setup_priors_extended(opts)
% Build DCM-style priors for the extended TCM.
% opts.free_mode: 'per-receptor' (default) | 'per-edge-fixed-mix' | 'class-tied'
% opts.free_patterns: cellstr of regex on labels to free
% opts.var_w: SD for free weights (log-space), default 0.5
% opts.var_tau: SD for time-constants (log-space), default 0.2
% opts.include_intrinsic: true (default)
% opts.include_mods: true (default)

if nargin<1, opts=struct; end
opts = setdef(opts,'free_mode','per-receptor');
opts = setdef(opts,'free_patterns',{'^w_.*'});   % free all syn weights by default
opts = setdef(opts,'var_w',0.5);
opts = setdef(opts,'var_tau',0.2);
opts = setdef(opts,'include_intrinsic',true);
opts = setdef(opts,'include_mods',true);

% --- Registries ---
Reg = build_default_registries();
ns = 1;                       % number of sources (areas) you’ll use for compiling
Comp = compile_tcm(Reg, ns);  %#ok<NASGU> % (compiled maps if needed later)

% --- Collect labels ---
labels = {};
% (i) Synaptic edge weights (per receptor)
labels = [labels; Reg.Syn.param];  % like 'w_SP_SI_AMPA'

% (ii) Time constants (taus): KE, KN, KI, KB, KM, KH, KKA, KSK, KNaP
taus = fieldnames(Reg.Params.tau);
labels = [labels; taus];

% (iii) Background/exogenous baseline (BE)
labels = [labels; {'BE'}];

% (iv) Intrinsic gains (optional)
if opts.include_intrinsic
    % e.g., leak/intrinsic channel gains per pop if you expose them
    % Add your own (gL_*, gM_*, gHCN_* ...) if you keep them parametric
    % For now we keep these minimal and let channels be driven by drive/gain postsyn mods
end

% (v) Modulators (optional): names in Reg.Mods
if opts.include_mods
    labels = [labels; {Reg.Mods.name}'];
end

labels = labels(:);
nP = numel(labels);

% --- Prior means (log-space)
pE = zeros(nP,1);
% Synaptic weights neutral:
pE(startsWith(labels,'w_')) = log(1.0);
% Taus: initialise to log of registry defaults
for k = 1:numel(taus)
    idx = strcmp(labels, taus{k});
    pE(idx) = log(Reg.Params.tau.(taus{k}));
end
% Baseline BE
pE(strcmp(labels,'BE')) = log(0.8);
% Modulator sensitivities (small)
if opts.include_mods
    for m = 1:numel(Reg.Mods)
        pE(strcmp(labels, Reg.Mods(m).name)) = log(0.1);
    end
end

% --- Prior covariance (fix-all then free selection)
pC = zeros(nP,nP);

free_mask = false(nP,1);
switch lower(opts.free_mode)
 case 'per-receptor'
   for r = 1:numel(opts.free_patterns)
     free_mask = free_mask | ~cellfun('isempty', regexp(labels, opts.free_patterns{r}, 'once'));
   end
 case 'per-edge-fixed-mix'
   % Free one scalar per (pre,post) and enforce AMPA/NMDA/KA or GABAA/B tying in your model
   % (We mark all matching receptors free here; the tying is enforced inside f/likelihood)
   for r = 1:numel(opts.free_patterns)
     free_mask = free_mask | ~cellfun('isempty', regexp(labels, opts.free_patterns{r}, 'once'));
   end
 case 'class-tied'
   % Free class-shared params by marking all within-class and tying in the model
   free_mask = ~cellfun('isempty', regexp(labels,'^w_.*_(AMPA|NMDA|KA)$','once')) | ...
               ~cellfun('isempty', regexp(labels,'^w_.*_(GABAA|GABAB)$','once'));
 otherwise
   error('Unknown free_mode: %s', opts.free_mode);
end

% Variances
pC(free_mask,free_mask) = diag( (opts.var_w.^2) * ones(nnz(free_mask),1) );

% Let taus move a bit (log-space) if desired:
tau_mask = ismember(labels, taus);
pC(tau_mask, tau_mask) = pC(tau_mask, tau_mask) + diag( (opts.var_tau.^2) * ones(nnz(tau_mask),1) );

% Small variance for mods (if included)
if opts.include_mods
    mod_mask = false(nP,1);
    for m = 1:numel(Reg.Mods)
        mod_mask = mod_mask | strcmp(labels,Reg.Mods(m).name);
    end
    pC(mod_mask, mod_mask) = pC(mod_mask, mod_mask) + diag( (0.2.^2) * ones(nnz(mod_mask),1) );
end

% --- Output
Pri.labels = labels;
Pri.pE     = pE;
Pri.pC     = pC;

% Convenience mappers (vector<->struct) matching label names
Pri.vec_to_struct = @(v) vec_to_struct(v, labels);
Pri.struct_to_vec = @(S) struct_to_vec(S, labels);

end

% ===== helpers =====
function s = setdef(s,f,v), if ~isfield(s,f)||isempty(s.(f)), s.(f)=v; end, end

function S = vec_to_struct(v, labels)
S = struct();
for k=1:numel(labels)
    S.(matlab.lang.makeValidName(labels{k})) = v(k);
end
end

function v = struct_to_vec(S, labels)
v = zeros(numel(labels),1);
for k=1:numel(labels)
    nm = matlab.lang.makeValidName(labels{k});
    if isfield(S,nm), v(k) = S.(nm); end
end
end
