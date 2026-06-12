function [M, pE, pC, x0] = make_tc_ngnmm_model(Ns)
% make_tc_ngnmm_model
% -------------------------------------------------------------------------
% Build an 8-population thalamo-cortical next-generation neural mass model
% in a DCM-compatible state format.
%
% This is the model as described in:
% Coombes S (2023) Next generation neural population models. 
% Front. Appl. Math. Stat. 9:1128224. doi: 10.3389/fams.2023.1128224
%
% and following the architecture of our previous 8-population
% thalamo-cortical model:
% Berndt, Singh & Shaw (2025) Restoring Synaptic Balance in Schizophrenia: 
% Insights From a Thalamo-Cortical Conductance-Based Model, 
% Schizophrenia Bulletin, 2025;, sbaf149, https://doi.org/10.1093/schbul/sbaf149
%
% This model implements a QIF/theta-neuron mean-field style neural mass
% system with alpha-function synaptic conductances.  It is designed to be
% used with DCM-style spectral fitting routines, where the model state is
% linearised around an operating point and passed to a Laplace-domain
% transfer function.
%
% Usage
% -------------------------------------------------------------------------
%   [M, pE, pC, x0] = make_tc_ngnmm_model(Ns)
%
% Inputs
% -------------------------------------------------------------------------
%   Ns      number of sources/regions. Default = 1.
%
% Outputs
% -------------------------------------------------------------------------
%   M       model structure containing dimensions, population labels,
%           connectivity masks, reversal potentials, state template and
%           model function handle.
%
%   pE      prior expectation structure. This is typically assigned to
%           DCM.M.pE.
%
%   pC      prior covariance structure. This is typically assigned to
%           DCM.M.pC.
%
%   x0      initial state array, size Ns x Npop x Nstate.
%
%
% Model populations
% -------------------------------------------------------------------------
% The model contains 8 populations:
%
%   1 SS    superficial stellate
%   2 SP    superficial pyramidal
%   3 SI    superficial inhibitory
%   4 DP    deep pyramidal
%   5 DI    deep inhibitory
%   6 TP    thalamic projection / thalamic pyramidal-like
%   7 RT    reticular thalamic inhibitory
%   8 RL    relay thalamic
%
% Population type is stored in:
%
%   M.popType
%
% with +1 for excitatory/projection-like populations and -1 for inhibitory
% populations.
%
%
% State variables
% -------------------------------------------------------------------------
% The state array has size:
%
%   x0 = Ns x Npop x Nstate
%
% where Nstate = 4:
%
%   x(:,:,1) = R    population firing rate
%   x(:,:,2) = V    mean membrane voltage / QIF voltage variable
%   x(:,:,3) = g    lumped incoming synaptic conductance
%   x(:,:,4) = h    conductance derivative, giving alpha-synapse dynamics
%
% The model function returns vectorised derivatives in the same ordering as
% x(:), making it compatible with SPM/DCM-style vectorisation.
%
%
% Intrinsic connectivity
% -------------------------------------------------------------------------
% Intrinsic connectivity is encoded by:
%
%   M.Kmask
%
% This is an Npop x Npop binary matrix. Rows are target populations and
% columns are source/pre-synaptic populations:
%
%   M.Kmask(target, source)
%
% A value of 1 means that the source population contributes to the target
% population's incoming conductance drive.  A value of 0 fixes that edge to
% a near-zero value by setting the prior covariance of the corresponding
% logKappa entry to zero.
%
%
% Parameters in DCM.M.pE / pE
% -------------------------------------------------------------------------
%
% pE.logDelta
%   1 x Npop vector.
%   Log half-width / dispersion parameter of the QIF population.
%   In the dynamics this becomes:
%
%       Delta = exp(P.logDelta)
%
%   Larger Delta increases the baseline drive to the firing-rate equation
%   and can make populations more excitable or broaden their response.
%
%
% pE.eta
%   1 x Npop vector.
%   Baseline excitability / tonic drive of each population in the QIF
%   voltage equation.
%
%   More positive eta pushes a population closer to active firing regimes.
%   More negative eta makes the population more subthreshold or quiescent.
%   In this implementation, excitatory/projection populations are given a
%   less negative default prior than inhibitory populations.
%
%
% pE.logAlpha
%   1 x Npop vector.
%   Log alpha-synapse rate for each target population:
%
%       alpha = exp(P.logAlpha)
%
%   This controls the speed of the lumped conductance kernel for each
%   target population:
%
%       dg/dt = h
%       dh/dt = alpha^2 * (drive_g - g) - 2*alpha*h
%
%   Larger alpha gives faster synaptic/conductance dynamics and can support
%   higher-frequency modes.  Smaller alpha gives slower conductance dynamics
%   and favours lower-frequency modes.
%
%
% pE.logKappa
%   Npop x Npop matrix.
%   Log intrinsic coupling strength:
%
%       Kappa = exp(P.logKappa)
%
%   Rows are target populations and columns are source populations:
%
%       Kappa(target, source)
%
%   The conductance drive to each target population is:
%
%       drive_g(target) = sum_source Kappa(target, source) * R(source)
%
%   Entries outside M.Kmask are set close to zero and are usually fixed by
%   pC.logKappa = 0 for absent edges.
%
%   Increasing an excitatory source-to-target Kappa usually strengthens
%   AMPA-like depolarising drive. Increasing an inhibitory source-to-target
%   Kappa usually strengthens GABA-like hyperpolarising drive, because the
%   reversal potential is determined by the source population type.
%
%
% pE.Erev
%   Npop x Npop matrix.
%   Reversal potential for each target x source connection.
%
%   Rows are target populations and columns are source populations:
%
%       Erev(target, source)
%
%   By default, excitatory/projection sources have Erev = +60 and
%   inhibitory sources have Erev = -90.  The model computes a
%   conductance-weighted effective reversal potential for each target
%   population:
%
%       Eeff(target) =
%           sum_source Kappa(target,source)*Erev(target,source)
%           ----------------------------------------------------
%                    sum_source Kappa(target,source)
%
%   Reversal potentials are usually fixed initially by setting pC.Erev = 0.
%
%
% pE.logQIFRate
%   Scalar.
%   Log global timescale/rate multiplier for the QIF firing-rate and voltage
%   equations:
%
%       qifRate = exp(P.logQIFRate)
%
%   This scales the R and V dynamics. Larger values make the QIF membrane
%   and rate dynamics faster, which can shift resonant modes upward in
%   frequency. Smaller values slow the intrinsic QIF dynamics.
%
%
% pE.inputGain
%   1 x Npop vector.
%   External input gain per population.
%
%   If a non-empty input u is supplied to the model function, the first
%   element of u is multiplied by inputGain and added to the voltage equation
%   of each population.  For resting CSD fits this often plays a minor role,
%   but it keeps the model compatible with DCM-style exogenous inputs.
%
%
% pE.d
%   1 x 3 vector.
%   Legacy delay parameter retained for compatibility with older TCM /
%   Laplace-transfer-function wrappers.  The first and third entries are
%   log-delay-like terms in the older convention.
%
%
% pE.logDelay
%   Scalar.
%   Log base delay in seconds:
%
%       baseDelay = exp(P.logDelay)
%
%   In tc_ngnmm_fx, this delay is placed on the path from pre-synaptic
%   firing rate R(source) to the post-synaptic conductance acceleration
%   h(target).  In other words, delays affect the R_pre -> h_post
%   dependencies that generate delayed synaptic conductance drive.
%
%
% pE.L
%   1 x Npop vector.
%   Optional legacy observer/source-leadfield-style field retained for older
%   wrappers that check for its existence.  In this implementation it is not
%   the main observation model; observation weights are encoded by pE.J.
%
%
% pE.J
%   numel(x0) x 1 vector.
%   Log observation weight for each vectorised state in x0(:).
%
%   There is one J parameter per state variable.  The model therefore allows
%   different hidden states and populations to contribute differently to the
%   observed spectrum.  By default, voltage states have weak visibility, with
%   stronger prior visibility on pyramidal and relay populations. Some
%   conductance states are also given weak prior visibility.
%
%   The corresponding linear observation weight is typically exp(P.J),
%   depending on the transfer-function / observation code.
%
%
% Parameters in DCM.M.pC / pC
% -------------------------------------------------------------------------
% pC has the same field structure as pE and controls prior uncertainty.
%
%   pC.field = 0      fixes that parameter at its prior expectation.
%   pC.field > 0      allows that parameter to vary during fitting.
%
% Key defaults:
%
%   pC.logDelta     permits population dispersion to vary.
%   pC.eta          permits baseline excitability to vary.
%   pC.logAlpha     permits synaptic timescales to vary.
%   pC.logKappa     permits only M.Kmask-present edges to vary.
%   pC.Erev         fixes reversal potentials initially.
%   pC.inputGain    usually weak or fixed, except selected populations.
%   pC.logQIFRate   permits the global QIF timescale to vary.
%   pC.logDelay     permits delay to vary weakly.
%   pC.J            controls which hidden states are observable.
%
%
% Frequency/mode interpretation
% -------------------------------------------------------------------------
% The main parameters controlling the linearised spectral modes are:
%
%   logAlpha       population-specific synaptic/conductance speed
%   logKappa       intrinsic coupling strengths
%   eta            baseline excitability / operating point
%   logDelta       population dispersion / excitability
%   logQIFRate     global QIF R/V timescale
%   logDelay       delayed synaptic drive
%
% In practice:
%
%   Increasing logAlpha or logQIFRate tends to make dynamics faster and can
%   help produce higher-frequency modes.
%
%   Changing logKappa alters loop gain and damping. Stronger recurrent
%   excitation can sharpen or destabilise modes, while stronger inhibition
%   can shift, damp, or stabilise modes depending on the loop.
%
%   Changing eta and logDelta moves the fixed point, which can substantially
%   change the eigenmodes of the linearised system.
%
%
% Notes
% -------------------------------------------------------------------------
% This model is deliberately pragmatic and DCM-friendly rather than a fully
% biophysical conductance model. It is intended to provide a compact
% next-generation neural mass approximation to an 8-population
% thalamo-cortical architecture suitable for spectral fitting.
%
% AS2026

if nargin < 1 || isempty(Ns)
    Ns = 1;
end

% ---------------------------------------------------------------------
% Dimensions and population labels
% ---------------------------------------------------------------------
M.Ns      = Ns;
M.Npop    = 8;
M.Nstates = 4;
M.pop     = {'SS','SP','SI','DP','DI','TP','RT','RL'};

Np = M.Npop;
Nx = M.Nstates;

% Population type: +1 excitatory/projection, -1 inhibitory
% SS, SP, DP, TP, RL treated as excitatory/projection-like;
% SI, DI, RT treated as inhibitory.
M.popType = [1 1 -1 1 -1 1 -1 1];

% ---------------------------------------------------------------------
% Intrinsic TCM-style target x source connectivity mask
% ---------------------------------------------------------------------
% Rows = target population, columns = source/pre-synaptic population.
% This is a pragmatic NGNMM version of the 8-pop TCM architecture.
% It can be tightened/altered later to match your exact tc_hilge2 graph.
Kmask = zeros(Np,Np);

SS = 1; SP = 2; SI = 3; DP = 4; DI = 5; TP = 6; RT = 7; RL = 8;

% Local cortical superficial circuit
Kmask(SS,SS) = 1;    % SS -> SS
Kmask(SP,SS) = 1;    % SS -> SP
Kmask(SI,SS) = 1;    % SS -> SI
Kmask(SS,SI) = 1;    % SI -> SS
Kmask(SP,SI) = 1;    % SI -> SP
Kmask(SI,SI) = 1;    % SI -> SI
Kmask(SI,SP) = 1;    % SP -> SI
Kmask(SP,SP) = 1;    % SP -> SI

% Superficial/deep cortical interactions
Kmask(DP,SP) = 1;    % SP -> DP
Kmask(DI,SP) = 1;    % SP -> DI
Kmask(SP,DP) = 1;    % DP -> SP
Kmask(DI,DP) = 1;    % DP -> DI
Kmask(DP,DI) = 1;    % DI -> DP
Kmask(DI,DI) = 1;    % DI -> DI

% Cortico-thalamic and thalamo-cortical loops
Kmask(TP,DP) = 1;    % DP -> TP
Kmask(RT,DP) = 1;    % DP -> RT
Kmask(RL,DP) = 1;    % DP -> RL
Kmask(SS,RL) = 1;    % RL -> SS
Kmask(SP,RL) = 1;    % RL -> SP
Kmask(SI,RL) = 1;    % RL -> SI

% Thalamic microcircuit
Kmask(RT,TP) = 1;    % TP -> RT
Kmask(RL,TP) = 1;    % TP -> RL
Kmask(TP,RL) = 1;    % RL -> TP
Kmask(RT,RL) = 1;    % RL -> RT
Kmask(TP,RT) = 1;    % RT -> TP
Kmask(RL,RT) = 1;    % RT -> RL
Kmask(RT,RT) = 1;    % RT -> RT

M.Kmask = Kmask;

% Reversal potential matrix, target x source.  The reversal is assigned
% by source population type: excitatory sources AMPA-like, inhibitory
% sources GABA-like.
Erev = 60 * ones(Np,Np);
Erev(:, M.popType < 0) = -90;
M.Erev = Erev;

% Population-specific default alpha rates, in Hz-like units.
% These are deliberately conservative starting values.
alpha0 = [100 100 50 100 50 80 40 80];
alpha0 = [60 60 40 60 40 50 35 50];
alpha0 = [25 30 20 25 20 25 18 25];
alpha0 = [45 55 35 45 35 45 30 45];
M.alpha0 = alpha0;

% % ---------------------------------------------------------------------
% % Initial states: Ns x Npop x Nstate
% % ---------------------------------------------------------------------
x0 = zeros(Ns,Np,Nx);
x0(:,:,1) = 0.10;     % R
x0(:,:,2) = -1.00;    % V
x0(:,:,3) = 0.01;     % g
x0(:,:,4) = 0.00;     % h


% ---------------------------------------------------------------------
% Priors
% ---------------------------------------------------------------------
pE = struct();
pC = struct();

% Heterogeneity / dispersion
pE.logDelta = log(0.10) * ones(1,Np);
pE.logDelta = log(0.20) * ones(1,Np);

% More subthreshold baseline drive
pE.eta = -3.0 * ones(1,Np);
pE.eta = -2.0 * ones(1,Np);

% Slightly less negative for relay/projection populations
%pE.eta([TP RL]) = -2.5;

pE.eta([SS SP DP TP RL]) = -1.5;
pE.eta([SI DI RT])       = -2.2;

% Slightly bias excitatory/projection populations toward activity and
% inhibitory populations toward restraint; keep broad enough for fitting.
pE.eta(M.popType > 0) = 0.0;
pE.eta(M.popType < 0) = -0.5;

pE.logAlpha = log(alpha0);
%pE.logAlpha = log(alpha0);

% Intrinsic target x source connectivity.  Absent edges are very small;
% present edges start moderately weak.  pC.logKappa below fixes absent
% edges by setting variance to zero.
pE.logKappa = log(1e-6 * ones(Np,Np));
pE.logKappa(Kmask == 1) = log(0.06);

% Slightly stronger canonical TCM edges
pE.logKappa(SP,SS) = log(0.10);
pE.logKappa(SI,SS) = log(0.08);
pE.logKappa(DP,SP) = log(0.08);
pE.logKappa(DI,DP) = log(0.08);
pE.logKappa(SS,RL) = log(0.08);
pE.logKappa(TP,RT) = log(0.08);
pE.logKappa(RL,RT) = log(0.08);

% Superficial cortical loop
pE.logKappa(SI,SS) = log(0.12);   % SS -> SI
pE.logKappa(SP,SI) = log(0.12);   % SI -> SP
pE.logKappa(SI,SP) = log(0.12);   % SI -> SP
pE.logKappa(SI,SI) = log(0.16);   % SI -> SP


% Deep cortical loop
pE.logKappa(DI,DP) = log(0.12);   % DP -> DI
pE.logKappa(DP,DI) = log(0.12);   % DI -> DP

% Thalamic loop
pE.logKappa(RL,RT) = log(0.12);   % RT -> RL
pE.logKappa(TP,RT) = log(0.10);   % RT -> TP
pE.logKappa(RT,TP) = log(0.08);   % TP -> RT, if allowed in Kmask

pE.Erev = Erev;

% QIF / Montbrio mean-field timescale, in s^-1
% 100 gives roughly two orders of magnitude faster R/V dynamics.
pE.logQIFRate = log(120); % 60, 80, 100, 120
pC.logQIFRate = 1/4;

% External input gain per population.  For CSD resting fits this will
% often not dominate, but it keeps the model DCM-compatible.
pE.inputGain = zeros(1,Np);
pE.inputGain(SS) = 1;
pE.inputGain(SP) = 0.5;

% Delay priors.  Existing Alex_LaplaceTFwDNew code often expects P.d.
% P.d is kept as a 3-vector for compatibility with older TCM code.
pE.d = [log(0.002) 0 log(0.002)];
pE.logDelay = log(0.002);

% Optional L field retained for older wrappers that check for it.
pE.L = zeros(1,Np);
pE.L(SP) = 1;

% Prior covariances
pC.logDelta = 1/4  * ones(1,Np);
pC.eta      = 1/2  * ones(1,Np);
pC.logAlpha = 1/4  * ones(1,Np);

pC.logKappa = zeros(Np,Np);
pC.logKappa(Kmask == 1) = 1/4;

% Usually keep reversals fixed initially.
pC.Erev = zeros(Np,Np);

pC.inputGain = zeros(1,Np);
pC.inputGain([RT RL]) = 1/8;

pC.d = [1 0 1] ./ 8;
pC.logDelay = 1/16;
pC.L = zeros(1,Np);

% ---------------------------------------------------------------------
% Observer weights: one per vectorised state in M.x(:)
% ---------------------------------------------------------------------
nx = numel(x0);

pE.J = -32 * ones(nx,1);
pC.J = zeros(nx,1);

for s = 1:Ns
    for p = 1:Np
        iV = sub2ind(size(x0),s,p,2);

        % Default weak voltage visibility
        pE.J(iV) = log(0.1);
        pC.J(iV) = 1/8;
    end

    % Stronger prior visibility on pyramidal / relay populations
    for p = [SP DP RL]
        iV = sub2ind(size(x0),s,p,2);

        pE.J(iV) = log(1.0);
        pC.J(iV) = 1;
    end
end

for s = 1:Ns
    for p = [SS SP DP RL]
        iG = sub2ind(size(x0),s,p,3);

        pE.J(iG) = log(0.05);
        pC.J(iG) = 1/8;
    end
end

% ---------------------------------------------------------------------
% Initial states
% ---------------------------------------------------------------------
x0 = zeros(Ns,Np,Nx);

% Start in a low-rate regime
R0 = 0.05 * ones(Ns,Np);

% Intrinsic drive implied by R0
Kappa0 = exp(pE.logKappa);
drive_g0 = zeros(Ns,Np);

for s = 1:Ns
    drive_g0(s,:) = (Kappa0 * R0(s,:)')';
end

% Make g initially self-consistent so dh is near zero
g0 = drive_g0;
h0 = zeros(Ns,Np);

% Choose V to approximately satisfy dR = 0:
% dR = Delta/pi + 2RV - gR = 0
Delta0 = exp(pE.logDelta(:))';
V0 = zeros(Ns,Np);

for s = 1:Ns
    for p = 1:Np
        V0(s,p) = 0.5 * (g0(s,p) - Delta0(p) / (pi * R0(s,p)));
    end
end

x0(:,:,1) = R0;
x0(:,:,2) = V0;
x0(:,:,3) = g0;
x0(:,:,4) = h0;



% ---------------------------------------------------------------------
% Attach to M
% ---------------------------------------------------------------------
M.x = x0;

% Use package-qualified handle if this file lives in +atcm.
try
    M.f = @atcm.tc_ngnmm_fx;
catch
    M.f = @tc_ngnmm_fx;
end

end
