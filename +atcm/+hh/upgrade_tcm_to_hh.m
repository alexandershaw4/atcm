function [M,P] = upgrade_tcm_to_hh(M,P,varargin)
%UPGRADE_TCM_TO_HH Expand tc_hilge2 model/prior structures for tc_hh_tcm.
%
%   [M,P] = upgrade_tcm_to_hh(M,P)
%   [M,P] = upgrade_tcm_to_hh(M,P,'statePriorVariance',1/16,'hhParamVariance',1/16)
%
% This upgrades an existing thalamo-cortical model from the old 7-state
% layout used by tc_hilge2 to the 10-state HH layout used by tc_hh_tcm.
%
% Old state layout, per source/population:
%   1 V, 2 AMPA, 3 GABA-A, 4 NMDA, 5 GABA-B, 6 M conductance, 7 H conductance
%
% New state layout, per source/population:
%   1 V, 2 AMPA, 3 GABA-A, 4 NMDA, 5 GABA-B,
%   6 pM, 7 rH, 8 mNa, 9 hNa, 10 nK
%
% The function also augments state priors and prior variances:
%   P.J       is expanded/created to length numel(M.x), with new HH states
%             initialised at zero in log/deviation space.
%   M.pC.J    is expanded/created to the same length, with configurable prior
%             variance for the added HH gates.
%
% It also adds conservative HH-specific parameter priors to P, with matching
% M.pC entries when absent.
%
% Typical use:
%   [M,P] = upgrade_tcm_to_hh(M,P);
%   M.f   = @tc_hh_tcm;
%   assert(numel(P.J) == numel(M.x));
%   assert(numel(M.pC.J) == numel(M.x));
%
% This function preserves wrapped parameter structures of the form P.p.

% -------------------------------------------------------------------------
% Options
% -------------------------------------------------------------------------
statePriorVariance = 1/16;   % prior variance for newly added HH state ICs
hhParamVariance    = 1/16;   % prior variance for newly added HH parameters
oldStateVariance   = 0;      % used only if M.pC.J is absent entirely

if mod(numel(varargin),2) ~= 0
    error('Optional arguments must be supplied as name/value pairs.');
end
for k = 1:2:numel(varargin)
    name = lower(varargin{k});
    val  = varargin{k+1};
    switch name
        case 'statepriorvariance'
            statePriorVariance = val;
        case 'hhparamvariance'
            hhParamVariance = val;
        case 'oldstatevariance'
            oldStateVariance = val;
        otherwise
            error('Unknown option: %s',varargin{k});
    end
end

if nargin < 2 || isempty(P)
    P = struct();
end

wrappedP = isstruct(P) && isfield(P,'p') && isstruct(P.p);
if wrappedP
    Pcore = P.p;
else
    Pcore = P;
end

% -------------------------------------------------------------------------
% Dimensions and state expansion
% -------------------------------------------------------------------------
if ~isfield(M,'x') || isempty(M.x)
    error('M.x is required.');
end

ns  = size(M.x,1);
np  = size(M.x,2);
nk0 = size(M.x,3);

if np ~= 8
    error('Expected M.x to have 8 populations, but size(M.x,2) is %d.',np);
end

oldNumel = numel(M.x);
oldNkForPriors = nk0;

if nk0 < 10
    xold = M.x;
    xnew = zeros(ns,np,10);

    % Preserve voltage and ligand-gated synaptic conductance states.
    ncopy = min(nk0,5);
    if ncopy > 0
        xnew(:,:,1:ncopy) = xold(:,:,1:ncopy);
    end

    % Use a physiologically sensible starting voltage if the original model
    % has zero initial voltages.
    V = xnew(:,:,1);
    if all(abs(V(:)) < 1e-12)
        V(:) = -65;
        xnew(:,:,1) = V;
    end

    % Initialise M/H gates from their steady-state activation curves, rather
    % than copying old M/H conductance states into activation variables.
    VhalfM = local_expand_linear(local_field(Pcore,'VhalfM',-35),ns,np);
    kM     = local_expand_linear(local_field(Pcore,'kM',10),ns,np);
    pMinf  = 1 ./ (1 + exp(-(V - VhalfM) ./ max(abs(kM),1e-6)));

    VhalfH = local_expand_linear(local_field(Pcore,'VhalfH',-75),ns,np);
    kH     = local_expand_linear(local_field(Pcore,'kH',-7),ns,np);
    rHinf  = 1 ./ (1 + exp((V - VhalfH) ./ max(abs(kH),1e-6)));

    Vsh = local_expand_linear(local_field(Pcore,'VshiftHH',0),ns,np);
    [am,bm,ah,bh,an,bn] = local_hh_rates(V + Vsh);

    xnew(:,:,6)  = local_clip01(pMinf);
    xnew(:,:,7)  = local_clip01(rHinf);
    xnew(:,:,8)  = local_clip01(am ./ max(am + bm,eps));
    xnew(:,:,9)  = local_clip01(ah ./ max(ah + bh,eps));
    xnew(:,:,10) = local_clip01(an ./ max(an + bn,eps));

    M.x = xnew;
else
    % Already upgraded. Keep all existing states, including any states beyond
    % the first 10, but still patch priors below.
    V = M.x(:,:,1);
end

newNumel = numel(M.x);
nkNew    = size(M.x,3);

% -------------------------------------------------------------------------
% Expand P.J, the log/deviation-space initial-condition parameter vector
% -------------------------------------------------------------------------
if isfield(Pcore,'J') && ~isempty(Pcore.J)
    Pcore.J = local_expand_state_vector(Pcore.J,ns,np,nkNew,0);
else
    Pcore.J = zeros(newNumel,1);
end

% -------------------------------------------------------------------------
% Expand M.pC.J, the prior variance for P.J
% -------------------------------------------------------------------------
if ~isfield(M,'pC') || isempty(M.pC)
    M.pC = struct();
end

if isstruct(M.pC)
    if isfield(M.pC,'J') && ~isempty(M.pC.J)
        M.pC.J = local_expand_state_vector(M.pC.J,ns,np,nkNew,statePriorVariance);
    else
        Jvar = oldStateVariance * ones(ns,np,nkNew);
        if nkNew >= 6
            Jvar(:,:,6:nkNew) = statePriorVariance;
        end
        M.pC.J = Jvar(:);
    end
else
    % If M.pC was numeric, do not destroy it. Add a separate M.pC_J fallback
    % so downstream code can still find the state prior variances.
    if isfield(M,'pC_J') && ~isempty(M.pC_J)
        M.pC_J = local_expand_state_vector(M.pC_J,ns,np,nkNew,statePriorVariance);
    else
        Jvar = oldStateVariance * ones(ns,np,nkNew);
        if nkNew >= 6
            Jvar(:,:,6:nkNew) = statePriorVariance;
        end
        M.pC_J = Jvar(:);
    end
end

% -------------------------------------------------------------------------
% Add conservative HH-specific parameter priors and matching variances
% -------------------------------------------------------------------------
% Log-scaled maximal conductance and rate parameters used by tc_hh_tcm.
[Pcore,M] = local_add_prior(Pcore,M,'gNa',      log([8 10 4 10 4 8 5 8]),       hhParamVariance);
[Pcore,M] = local_add_prior(Pcore,M,'gK',       log([3 4 2 4 2 4 3 4]),          hhParamVariance);
[Pcore,M] = local_add_prior(Pcore,M,'gL',       log(ones(1,8)),                  hhParamVariance/4);
[Pcore,M] = local_add_prior(Pcore,M,'gM',       log([0.4 0.6 0.3 0.8 0.3 1.4 0.3 1.4]), hhParamVariance);
[Pcore,M] = local_add_prior(Pcore,M,'gH',       log([0 0 0 0 0 1.2 0 1.6] + 1e-6), hhParamVariance);
[Pcore,M] = local_add_prior(Pcore,M,'phiHH',    log(0.15),                      hhParamVariance);

% Linear voltage/gating parameters used by tc_hh_tcm.
[Pcore,M] = local_add_prior(Pcore,M,'VshiftHH', 0,       4);
[Pcore,M] = local_add_prior(Pcore,M,'VhalfM',   -35,     9);
[Pcore,M] = local_add_prior(Pcore,M,'kM',       10,      4);
[Pcore,M] = local_add_prior(Pcore,M,'tauM',     log(1),  hhParamVariance);
[Pcore,M] = local_add_prior(Pcore,M,'VhalfH',   -75,     9);
[Pcore,M] = local_add_prior(Pcore,M,'kH',       -7,      4);
[Pcore,M] = local_add_prior(Pcore,M,'tauH',     log(1),  hhParamVariance);

% Linear masks and optional reversal-potential fields. Keep these effectively
% fixed unless you deliberately open them by increasing M.pC.<field>.
[Pcore,M] = local_add_prior(Pcore,M,'NaPop',    [1 1 0.45 1 0.45 0.8 0.45 0.8], 0);
[Pcore,M] = local_add_prior(Pcore,M,'KPop',     [1 1 0.8  1 0.8  1   0.8  1],    0);
[Pcore,M] = local_add_prior(Pcore,M,'MPop',     [0.4 0.6 0.2 0.8 0.2 1.0 0.1 1.0], 0);
[Pcore,M] = local_add_prior(Pcore,M,'HPop',     [0 0 0 0 0 1 0 1],              0);
[Pcore,M] = local_add_prior(Pcore,M,'ENaHH',    55,      0);
[Pcore,M] = local_add_prior(Pcore,M,'EKHH',     -90,     0);
[Pcore,M] = local_add_prior(Pcore,M,'EHh',      -30,     0);

% Store a compact description of the new layout for later sanity checks.
M.hh_state_names = {'V','gA_AMPA','gG_GABAA','gN_NMDA','gB_GABAB', ...
                    'pM_Kv7','rH_HCN','mNa','hNa','nK'};
M.f = @atcm.hh.tc_hh_tcm;

% Re-wrap P if necessary.
if wrappedP
    P.p = Pcore;
else
    P = Pcore;
end

% Final sanity checks.
if wrappedP
    Jcheck = P.p.J;
else
    Jcheck = P.J;
end
if numel(Jcheck) ~= newNumel
    error('P.J expansion failed: numel(P.J)=%d but numel(M.x)=%d.',numel(Jcheck),newNumel);
end
if isstruct(M.pC) && isfield(M.pC,'J') && numel(M.pC.J) ~= newNumel
    error('M.pC.J expansion failed: numel(M.pC.J)=%d but numel(M.x)=%d.',numel(M.pC.J),newNumel);
end

end

% =========================================================================
% Local utilities
% =========================================================================
function val = local_field(P,field,default)
if isstruct(P) && isfield(P,field) && ~isempty(P.(field))
    val = P.(field);
else
    val = default;
end
end

function X = local_expand_linear(x,ns,np)
if isscalar(x)
    X = repmat(x,ns,np);
elseif isequal(size(x),[1 np])
    X = repmat(x,ns,1);
elseif isequal(size(x),[ns np])
    X = x;
elseif numel(x) == np
    X = repmat(reshape(x,1,np),ns,1);
elseif numel(x) == ns*np
    X = reshape(x,ns,np);
else
    error('Cannot expand parameter with size %s to ns x np.',mat2str(size(x)));
end
end

function vout = local_expand_state_vector(vin,ns,np,nkTarget,fillValue)
% Expand a state-sized vector/array from ns*np*nkOld to ns*np*nkTarget.
vin = full(vin);
if numel(vin) == ns*np*nkTarget
    vout = vin(:);
    return
end
if mod(numel(vin),ns*np) ~= 0
    error('State-sized parameter has %d elements, which is not divisible by ns*np=%d.',numel(vin),ns*np);
end
nkOld = numel(vin)/(ns*np);
Xold = reshape(vin,ns,np,nkOld);
Xnew = fillValue * ones(ns,np,nkTarget);
Xnew(:,:,1:min(nkOld,nkTarget)) = Xold(:,:,1:min(nkOld,nkTarget));
vout = Xnew(:);
end

function [P,M] = local_add_prior(P,M,field,defaultValue,varianceValue)
% Add P.<field> and M.pC.<field> only when absent. This keeps existing priors
% untouched but ensures the HH model has all fields it expects.
if ~isfield(P,field) || isempty(P.(field))
    P.(field) = defaultValue;
end
if ~isfield(M,'pC') || isempty(M.pC) || ~isstruct(M.pC)
    return
end
if ~isfield(M.pC,field) || isempty(M.pC.(field))
    M.pC.(field) = varianceValue * ones(size(P.(field)));
end
end

function [am,bm,ah,bh,an,bn] = local_hh_rates(V)
am = 0.1  .* local_vtrap(V + 40,10);
bm = 4.0  .* exp(-(V + 65)./18);
ah = 0.07 .* exp(-(V + 65)./20);
bh = 1.0  ./ (1 + exp(-(V + 35)./10));
an = 0.01 .* local_vtrap(V + 55,10);
bn = 0.125.* exp(-(V + 65)./80);
end

function y = local_vtrap(x,yden)
z = x ./ yden;
y = x ./ (1 - exp(-z));
small = abs(z) < 1e-6;
y(small) = yden .* (1 + z(small)./2);
end

function X = local_clip01(X)
X = min(max(X,0),1);
end
