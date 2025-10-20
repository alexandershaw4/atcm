function [f,J,Q,D] = tc_twocmp_stp(x,u,P,M)
% TC_TWOCMP_STP  Thalamo–cortical neural-mass with 2-compartment pyramids and presynaptic STP.
%
%   [f,J,Q,D] = atcm.tc_twocmp_stp(x,u,P,M)
%
% Implements a conductance-based canonical thalamo–cortical circuit with:
%   (i)  explicit PRE/POST split: presynaptic short-term plasticity (STP) per
%        presynaptic population (Tsodyks–Markram: resources R and use u);
%   (ii) two-compartment pyramids for superficial & deep pyramidal populations
%        (soma Vs, apical dendrite Vd) with axial coupling g_c and compartment-
%        specific synapse placement (E→dendrite, I→soma by default).
%
% Compatible with the original tc_hilge calling pattern and SPM-style DCM
% wrappers. Supports delay integration via the returned delay operator Q.
%
% -------------------------------------------------------------------------
% INPUTS
%   x : hidden states, sized [ns × 8 × nk].  Expected nk = 10 (see below).
%       If nk < 10, the function can be adapted to auto-pad (see code).
%   u : exogenous input(s). Scalar or vector; routed to RL (pop 8) and TP (6).
%   P : parameter struct. Standard fields from tc_hilge are respected:
%         A{1..5}    : extrinsic AMPA gains   (log space)
%         AN{1..5}   : extrinsic NMDA gains   (log space)
%         H(np×np×ns): intrinsic gains (AMPA/GABA template) (log space)
%         Hn(np×np×ns): intrinsic NMDA gains (log space)
%         T(np×6)    : channel time constants (log space), columns:
%                        [AMPA  GABAa  NMDA  GABAb  M  H]
%         C(ns×1)    : input gain(s) (log space)
%         CV(1×8)    : membrane capacitance scalers (log space)
%         pr         : firing midpoint shift (log space; VR = -52*exp(pr))
%         scale(4×1) : optional scaling of {GEa,GEn,GIa,GIb}
%         D(2×1)     : intrinsic delay scales for within-source (1) and
%                      across-source (2) interactions (log space)
%         D0(2×1)    : thalamo-cortical delay scales: [CT, TC] (log space)
%       Additional fields introduced here:
%         prel(8×1)  : baseline release probability per presyn pop (log)
%         tauR(8×1)  : STP resource recovery time constants (s)
%         tauU(8×1)  : STP facilitation time constants (s)
%         U0(8×1)    : baseline 'use' (logit space; model applies logistic)
%         gc         : axial coupling (log space; default ~3 mS)
%         w_dend(1×2): [AMPA_to_dend, NMDA_to_dend] in [0,1] (default [0.8 0.9])
%         w_soma(1×2): [GABAa_to_soma, GABAb_to_soma] in [0,1] (default [0.9 0.9])
%         IncludeMH  : 1/0 to include M/H currents on pops 6 & 8 (default 1)
%         E          : background excitatory drive (log space; default 0)
%
%   M : model struct with at least:
%         M.x : state template (used for ns,np,nk and Jacobian sizing)
%         (Optional fields used by your pipeline are passed through)
%
% -------------------------------------------------------------------------
% STATE LAYOUT (nk = 10 expected; per population)
%   1  Vs   : somatic membrane potential (mV)
%   2  gE   : AMPA conductance
%   3  gI   : GABA-A conductance
%   4  gN   : NMDA conductance (with Mg²⁺ block)
%   5  gB   : GABA-B conductance
%   6  gM   : M-current (optional; pops 6 & 8 by default)
%   7  gH   : H-current (optional; pops 6 & 8 by default)
%   8  Vd   : apical dendritic potential (used by pops 2 & 4; passive otherwise)
%   9  R    : STP resources (presynaptic, one per presyn population)
%   10 uSTP : STP 'use' (presynaptic)
%
% POPULATIONS (as tc_hilge):
%   1 SS (L4), 2 SP (L2/3), 3 SI (L2/3), 4 DP (L5), 5 DI (L5),
%   6 TP (L6), 7 RT (TRN), 8 RL (TC relay)
%
% -------------------------------------------------------------------------
% OUTPUTS
%   f : vectorised state derivatives, spm_vec-compatible with M.x
%   J : Jacobian dfdx (computed via SPM numerical differentiation if requested)
%   Q : delay operator,   Q = inv(I − D .* J)
%   D : state-space delay matrix (same size as J).  D combines:
%         • within-source, across-population delays (P.D(1))
%         • across-source delays (P.D(2))
%         • population-level thalamo-cortical delays:
%             CT (cortex→thal) and TC (thal→cortex) via P.D0
%       Units are seconds (negative sign convention as in tc_hilge).
%
% -------------------------------------------------------------------------
% MODEL NOTES
% • Presynaptic stage: x_pre = p_rel .* uSTP .* R .* m(Vs)
%   drives postsynaptic gates (AMPA/NMDA). Inhibition (GABAa/b) currently
%   follows firing rates m(Vs) (can be extended to presyn STP similarly).
% • Two-compartment pyramids (SP, DP): AMPA/NMDA primarily on dendrite,
%   GABAa/b on soma; set via w_dend and w_soma. Axial current ~ g_c(Vd−Vs).
% • Observation tip: for M/EEG, consider using a dipole-like readout
%   proportional to g_c*(Vd − Vs) from pyramidal populations.
%
% -------------------------------------------------------------------------
% DEFAULTS & RANGES (informal)
%   gc ≈ 3 mS (log-normal prior), w_dend ≈ [0.8 0.9], w_soma ≈ [0.9 0.9]
%   tauR ≈ 0.6 s, tauU ≈ 0.2 s, U0 ≈ 0.2, prel ≈ 0.6
%   T(:,*) follow tc_hilge scales: AMPA~2.2 ms, GABAa~5 ms, NMDA~100 ms,
%   GABAb~300 ms, M~160 ms, H~100 ms (all entered in log space).
%
% -------------------------------------------------------------------------
% EXAMPLES
%   % Switch to this model, expand state dims, and run once:
%   DCM.M.f = @atcm.tc_twocmp_stp;
%   if size(DCM.M.x,3) < 10
%       X = zeros(size(DCM.M.x,1), size(DCM.M.x,2), 10);
%       X(:,:,1:size(DCM.M.x,3)) = DCM.M.x;
%       X(:,:,8)  = X(:,:,1);    % Vd ~ Vs
%       X(:,:,9)  = 0.9;         % R
%       X(:,:,10) = 0.2;         % uSTP
%       DCM.M.x = X;
%   end
%   [f,J,Q,D] = DCM.M.f(DCM.M.x, 0, DCM.M.pE, DCM.M);
%
%   % Minimal prior additions (log/logit spaces):
%   DCM.M.pE.gc      = log(3);
%   DCM.M.pE.prel    = log(0.6)*ones(8,1);
%   DCM.M.pE.tauR    = 0.6*ones(8,1);
%   DCM.M.pE.tauU    = 0.2*ones(8,1);
%   DCM.M.pE.U0      = log(0.2/(1-0.2))*ones(8,1);   % logit
%   DCM.M.pE.w_dend  = [0.8 0.9];
%   DCM.M.pE.w_soma  = [0.9 0.9];
%
% -------------------------------------------------------------------------
% SEE ALSO
%   atcm.tc_hilge, atcm.Alex_LaplaceTFwD, spm_vec, spm_unvec, spm_diff
%
% ----------------------------------------------------------------------
% ADS2025

%gx = @(fld,def) (isfield(P,fld) && ~isempty(P.(fld))) * P.(fld) + (~isfield(P,fld) || isempty(P.(fld))) * def;
gx = @(fld,def) iff(isfield(P,fld) && ~isempty(P.(fld)), P.(fld), def);
function out = iff(cond, a, b), if cond, out = a; else, out = b; end, end

mg_switch = @(V) 1./(1 + 0.28*exp(-0.062*V));   % simple Mg2+ block
sig      = @(z) 1./(1+exp(-z));                 % logistic

% ----- basic dims & reshape ------------------------------------------------
ns = size(M.x,1);       % sources
np = size(M.x,2);       % populations (8)
nk = size(M.x,3);       % states per pop (expect 10 here)
x  = reshape(x,ns,np,nk);

% ----- indices (see header) -----------------------------------------------
iVs = 1; iGE = 2; iGI = 3; iGN = 4; iGB = 5; iGM = 6; iGH = 7; iVd = 8; iR = 9; iU = 10;

% ----- parameter defaults / transforms ------------------------------------
% Extrinsics
A  = cell(size(P.A)); AN = cell(size(P.A));
for k = 1:numel(P.A),  A{k}  = exp(P.A{k}); AN{k} = exp(P.AN{k}); end
C  = exp(gx('C',zeros(ns,1)));

% Intrinsics
G     = exp(full(gx('H' ,zeros(np,np,ns))));   % AMPA/GABA base (used with GEa/GIa)
Gn    = exp(full(gx('Hn',zeros(np,np,ns))));   % NMDA base (used with GEn)
scale = gx('scale',zeros(4,1));

% time constants
Tmat = gx('T', zeros(np,6));        % <-- np, not ns

KE = exp(-Tmat(:,1))*1000/2.2;        % AMPA
KI = exp(-Tmat(:,2))*1000/5;          % GABA-A
KN = exp(-Tmat(:,3))*1000/100;        % NMDA
KB = exp(-Tmat(:,4))*1000/300;        % GABA-B
KM = exp(-Tmat(:,5))*1000/160;        % M
KH = exp(-Tmat(:,6))*1000/100;        % H

% Voltages (mV)
VL = -70; VE =  60; VI = -90; VN = 10;  VB = -100;
VR = -52 * exp(gx('pr',0));                    % firing midpoint (allow shift)
CV = exp(gx('CV',zeros(1,8))).*[128*3 128 64 128 64 128 64 128*2]/1000; % per-pop C

% Optional M/H placement retained for compat
try; IncludeMH = gx('IncludeMH',1); catch; IncludeMH = 1; end
VM = -70; VH = -30;

% Connectivity templates from tc_hilge (lightly adapted)
[GEa,GEn,GIa,GIb,SA,SNMDA] = local_connectivity_templates(scale);

% ----- Two-compartment weights & axial coupling ---------------------------
gc      = exp(gx('gc',log(3)));             % axial coupling (mS) default ~3
wE_dend = min(max(gx('w_dend',[0.8 0.9]),0),1);   % [AMPA, NMDA] to dend
wI_soma = min(max(gx('w_soma',[0.9 0.9]),0),1);   % [GABAa, GABAb] to soma

% which pops are pyramids with two compartments
isPyr = false(1,8); isPyr([2 4]) = true;

% ----- Presynaptic STP (per PRESYN population) ----------------------------
% States R (resources), u (use). Per-population dynamics feed all its outputs.
U0   = sig(gx('U0',zeros(8,1)));     % baseline use
tauR = gx('tauR',0.6*ones(8,1));     % s
tauU = gx('tauU',0.2*ones(8,1));     % s
prel = exp(gx('prel',log(0.6)*ones(8,1)));  % baseline release prob per pop

% ----- firing nonlinearity -------------------------------------------------
FF = 1./(1 + exp(-exp(gx('S',zeros(1,8))).*(x(:,:,iVs)-VR)));
FF( x(:,:,iVs) >= VR ) = 1;          % clip above threshold region
FF( x(:,:,iVs) >= 30 ) = 0;          % refractory clip
m  = FF;                             % firing/output per population

% ----- extrinsics (as in original) ----------------------------------------
a       = zeros(ns,5); an = zeros(ns,5); 
a(:,1)  = A{1} * m(:,2);    an(:,1) = AN{1}*m(:,2);
a(:,2)  = A{2} * m(:,4);    an(:,2) = AN{2}*m(:,4);
a(:,3)  = A{3} * m(:,6);    an(:,3) = AN{3}*m(:,6);
a(:,4)  = A{4} * m(:,7);    an(:,4) = AN{4}*m(:,7);
a(:,5)  = A{5} * m(:,8);    an(:,5) = AN{5}*m(:,8);

BE = exp(gx('E',0))*0.8;             % background

% ----- allocate output container ------------------------------------------
f = zeros(ns,np,nk);

% ----- loop over sources ---------------------------------------------------
for is = 1:ns

    % ---- presyn STP dynamics (per presyn population) ----
    R  = x(is,:,iR);    % 1x8
    uS = x(is,:,iU);

    % presynaptic effective release x_pre for EACH presyn pop
    xpre = prel(:)'.*uS.*R.*m(is,:);     % 1x8

    % STP ODEs (Tsodyks–Markram)
    f(is,:,iR) = ((1 - R)./tauR' - uS.*R.*m(is,:));                     % dR/dt
    f(is,:,iU) = ((U0' - uS)./tauU' + U0'.*(1 - uS).*m(is,:));          % du/dt

    % ---- build post-synaptic drives via x_pre and templates ----
    % excitatory to post (AMPA/NMDA)
    E_AMPA = ((G(:,:,is).*GEa) * xpre(:))';      % 1x8
    E_NMDA = ((Gn(:,:,is).*GEn) * xpre(:))';     % 1x8
    % inhibitory to post (GABAa/GABAb) – use xpre from inhibitory presyn
    I_GABAA = ((G(:,:,is).*GIa) * m(is,:)' )';   % fast inhibit from rates
    I_GABAB = ((G(:,:,is).*GIb) * m(is,:)' )';   % slow inhibit from rates

    % Add extrinsic + background
    E_AMPA = (E_AMPA + BE + (SA   * a (is,:)')')' * 2;
    E_NMDA = (E_NMDA + BE + (SNMDA* an(is,:)')')' * 2;

    % Optionally route exogenous input to RL (8) and TP (6)
    dU = u(:)*C(is,1);
    if ~isempty(dU)
        if numel(dU) > 1
            E_AMPA(8) = E_AMPA(8) + dU(1);
            E_AMPA(6) = E_AMPA(6) + dU(2);
        else
            E_AMPA([8 1]) = E_AMPA([8 1]) + dU;
        end
    end

    % ---- conductance gating ODEs (post) ----
    % AMPA/NMDA driven by x_pre sums; GABAa/b by inhibitory m()
    f(is,:,iGE) = (E_AMPA' - x(is,:,iGE)) .* KE(is);
    f(is,:,iGN) = (E_NMDA' - x(is,:,iGN)) .* KN(is);
    f(is,:,iGI) = (I_GABAA - x(is,:,iGI)) .* KI(is);
    f(is,:,iGB) = (I_GABAB - x(is,:,iGB)) .* KB(is);

    % Optional intrinsic M/H
    if IncludeMH
        % simple mean-field currents proportional to firing for now
        Im = 4*[0 0 0 0 0 1 0 1] .* m(is,:);
        Ih = 4*[0 0 0 0 0 1 0 1] .* (1 - m(is,:));   % crude
        f(is,:,iGM) = (Im - x(is,:,iGM)) .* KM(is);
        f(is,:,iGH) = (Ih - x(is,:,iGH)) .* KH(is);
    end

    % ---- membrane equations -----------------------------------------------
    % SOMA currents
    GE = x(is,:,iGE); GI = x(is,:,iGI); GN = x(is,:,iGN); GB = x(is,:,iGB);
    GM = x(is,:,iGM); GH = x(is,:,iGH);
    Vs = x(is,:,iVs); Vd = x(is,:,iVd);

    % distribute synapses by compartment for pyramids (2 & 4)
    % soma: GABAa/b mostly; dend: AMPA/NMDA mostly
    GE_s = (1 - wE_dend(1)) * GE;   GE_d = wE_dend(1) * GE;
    GN_s = (1 - wE_dend(2)) * GN;   GN_d = wE_dend(2) * GN;
    GI_s = wI_soma(1) * GI;         GI_d = (1 - wI_soma(1)) * GI; % usually ~0 dend
    GB_s = wI_soma(2) * GB;         GB_d = (1 - wI_soma(2)) * GB; % usually ~0 dend

    % NMDA Mg2+ block on each compartment potential it acts on
    mgS = mg_switch(Vs);
    mgD = mg_switch(Vd);

    % Soma currents (per-pop)
    I_soma = (Vs - VL) .* (-1) ...                 % leak (GL=1)
           + GE_s .* (VE - Vs) ...
           + GI_s .* (VI - Vs) ...
           + GB_s .* (VB - Vs) ...
           + GN_s .* (VN - Vs) .* mgS ...
           + GM   .* (VM - Vs) ...
           + GH   .* (VH - Vs);

    % Dend currents
    I_dend = (Vd - VL) .* (-1) ...
           + GE_d .* (VE - Vd) ...
           + GI_d .* (VI - Vd) ...
           + GB_d .* (VB - Vd) ...
           + GN_d .* (VN - Vd) .* mgD;

    % Axial coupling: only for pyramids; others: set gc_eff=0
    gc_eff = zeros(1,8); gc_eff(isPyr) = gc;

    % dVs/dt and dVd/dt
    dVs = I_soma - gc_eff .* (Vs - Vd);
    dVd = I_dend - gc_eff .* (Vd - Vs);

    % For non-pyramids, collapse to single compartment (Vd passive to Vs)
    dVd(~isPyr) = 0;

    % divide by capacitance per-population
    f(is,:,iVs) = dVs ./ CV;
    f(is,:,iVd) = dVd ./ CV;

end

% ----- vectorise -----------------------------------------------------------
f = spm_vec(f);

% ----- Jacobian if requested ----------------------------------------------
if nargout < 2 || nargout == 5, J = []; return; end
J = spm_cat(spm_diff(M.f,spm_unvec(spm_vec(x),M.x),u,P,M,1));

% --- Delay operator (D) and Q = inv(I - D.*J) ---
if nargout >= 3
    % safe getters for delay params
    PD  = zeros(2,1);   if isfield(P,'D')  && ~isempty(P.D),  PD  = P.D(:);  end
    PD0 = zeros(2,1);   if isfield(P,'D0') && ~isempty(P.D0), PD0 = P.D0(:); end

    % Base intrinsic/self delays (seconds). Convention matches tc_hilge:
    % Dvec(1) applies within-source, across-population; Dvec(2) across sources.
    Dvec = [1 16];                          % ms base
    d    = - Dvec(:) .* exp(PD(1:numel(Dvec))) / 1000;   % -> seconds, negative for operator

    % Thalamo–cortical population delays (ms), with optional scaling in P.D0
    CTms = 8  * exp(PD0(1));   % cortex->thalamus
    TCms = 3  * exp(PD0(2));   % thalamus->cortex

    % population-level delay matrix (np x np), then expand to state-space
    Tc = zeros(np,np);
    Tc(1:6,[7 8]) = CTms;      % cortex -> thalamus
    Tc([7 8],1:6) = TCms;      % thalamus -> cortex
    Tc = -Tc / 1000;           % seconds, negative for operator

    % expand Tc to (ns*np*nk)-square
    Tc = kron(ones(nk,nk), kron(Tc, eye(ns)));

    % helpers to mark “same source” and “same population” in state indexing
    Sp = kron(ones(nk,nk), kron( eye(np), eye(ns)));  % same population & same source
    Ss = kron(ones(nk,nk), kron( ones(np), eye(ns))); % same source

    % Different-source vs same-source-different-population selectors
    Dp = ~Sp & ~Ss;   % different sources
    Ds = ~Sp &  Ss;   % same source, different population

    % Full delay matrix
    D  = d(2)*Dp + d(1)*Ds + Tc;

    % Delay operator
    N  = size(J,1);
    Q  = spm_inv(speye(N) - D.*J);
end

% -------------------------------------------------------------------------
% templates (mostly copied from your tc_hilge, lightly cleaned)
function [GEa,GEn,GIa,GIb,SA,SNMDA] = local_connectivity_templates(scale)
    GEa = [  0 0 0 0 0 2 0 2;
             4 4 0 0 0 0 0 0;
             0 2 0 0 0 0 0 0;
             0 2 0 0 0 0 0 0;
             0 0 0 2 0 0 0 0;
             0 0 0 8 0 0 0 0;
             0 0 0 0 0 0 0 2;
             2 0 0 0 0 2 0 0] * 2;
    GEn = GEa;
    si  = 8;
    GIa = [ si 0  8  0  0  0  0  0;
            0  si 64 0  0  0  0  0;
            0  0  si 0  0  0  0  0;
            0  0  0  12 4  0  0  0;
            0  0  16 0  4  0  0  0;
            0  0  0  0  8  4  0  0;
            0  0  0  0  0  0  12 0;
            0  0  0  0  0  0  4  si] * 2;
    GIb = GIa;
    if ~isempty(scale)
        GEa = GEa * exp(scale(1));
        GEn = GEn * exp(scale(2));
        GIa = GIa * exp(scale(3));
        GIb = GIb * exp(scale(4));
    end
    SA = [1 0 0 0 0;
          0 1 0 0 0;
          0 1 0 0 0;
          0 0 0 0 0;
          0 0 0 0 0;
          0 0 0 0 0;
          0 0 0 0 1;
          0 0 0 1 0]/8;
    SA(:,[3 4 5]) = 0;
    SNMDA = SA;
end
end
