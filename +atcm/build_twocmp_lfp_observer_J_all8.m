function [J,Wsrc,obs] = build_twocmp_lfp_observer_J_all8(M,P,weights,x0)
% BUILD_TWOCMP_LFP_OBSERVER_J_ALL8
% -------------------------------------------------------------------------
% Current-based LFP observation model for tc_twocmp_stp_all8.m.
%
% This replaces the older voltage/axial proxy
%       y ~ gc*(Vd - Vs)
% with a Mazzoni/Linden/Einevoll-style synaptic-current proxy:
%
%       D_p(t) = I_E,dend,p(t - tauE) - alpha * I_I,soma,p(t - tauI)
%
% and then mixes pyramidal/dipolar populations into a source signal:
%
%       y_source(t) = sum_p src_w(p) * D_p(t)
%
% The returned matrix J is the local linearisation d y / d vec(x), evaluated
% at x0 (default M.x). This is compatible with transfer-function / DCM code
% that expects a fixed observer matrix. Static J cannot itself implement a
% true temporal delay; the recommended delays are returned in obs.tauE and
% obs.tauI and can be applied in the frequency domain as exp(-s*tau).
%
% INPUTS
%   M.x      [ns x 8 x 10] state template / expansion point
%   P        parameter struct used by tc_twocmp_stp_all8
%   weights  optional struct:
%       .src          1x8 population source weights. Default SP=1, DP=0.6.
%       .alpha        inhibitory weight. Default 1.65.
%       .tauE         excitatory delay in seconds. Default 0.006.
%       .tauI         inhibitory delay in seconds. Default 0.
%       .useAMPA      include AMPA dendritic current. Default true.
%       .useNMDA      include NMDA dendritic current. Default true.
%       .useGABAA     include GABA-A soma current. Default true.
%       .useGABAB     include GABA-B soma current. Default true.
%       .normaliseRows normalise each source row by its L2 norm. Default false.
%       .projectLeadfield if true and M.L exists, return M.L*Wsrc. Default true.
%       .w_dend       override P.w_dend, rows=pops, cols=[AMPA NMDA].
%       .w_soma       override P.w_soma, rows=pops, cols=[GABAa GABAb].
%   x0       optional expansion point, same size as M.x
%
% OUTPUTS
%   J        nchan x numel(M.x) if M.L exists, otherwise ns x numel(M.x)
%   Wsrc     ns x numel(M.x) source-level observation Jacobian
%   obs      details, including current components at x0 and delay metadata
%
% STATE LAYOUT assumed from tc_twocmp_stp_all8:
%   1 Vs, 2 gE, 3 gI, 4 gN, 5 gB, 6 gM, 7 gH, 8 Vd, 9 R, 10 uSTP
%
% POPULATIONS:
%   1 SS, 2 SP, 3 SI, 4 DP, 5 DI, 6 TP, 7 RT, 8 RL/RC
% -------------------------------------------------------------------------

    if nargin < 3 || isempty(weights), weights = struct(); end
    if nargin < 4 || isempty(x0), x0 = M.x; end

    % ---- defaults --------------------------------------------------------
    if ~isfield(weights,'src') || isempty(weights.src)
        weights.src = zeros(1,8);
        weights.src(2) = 1.0;   % superficial pyramidal population
        weights.src(4) = 0.6;   % deep pyramidal population
    end
    if ~isfield(weights,'alpha')        || isempty(weights.alpha),        weights.alpha = 1.65;  end
    if ~isfield(weights,'tauE')         || isempty(weights.tauE),         weights.tauE = 0.006;  end
    if ~isfield(weights,'tauI')         || isempty(weights.tauI),         weights.tauI = 0.000;  end
    if ~isfield(weights,'useAMPA')      || isempty(weights.useAMPA),      weights.useAMPA = true; end
    if ~isfield(weights,'useNMDA')      || isempty(weights.useNMDA),      weights.useNMDA = true; end
    if ~isfield(weights,'useGABAA')     || isempty(weights.useGABAA),     weights.useGABAA = true; end
    if ~isfield(weights,'useGABAB')     || isempty(weights.useGABAB),     weights.useGABAB = true; end
    if ~isfield(weights,'normaliseRows')|| isempty(weights.normaliseRows),weights.normaliseRows = false; end
    if ~isfield(weights,'projectLeadfield') || isempty(weights.projectLeadfield), weights.projectLeadfield = true; end

    % ---- dimensions / indices -------------------------------------------
    ns = size(M.x,1); np = size(M.x,2); nk = size(M.x,3);
    if np ~= 8 || nk < 8
        error('Expected M.x to have size [ns x 8 x nk] with nk >= 8.');
    end
    x0 = reshape(x0,ns,np,nk);

    iVs = 1; iGE = 2; iGI = 3; iGN = 4; iGB = 5; iVd = 8;
    N   = numel(M.x);
    idx = @(is,ip,ik) sub2ind([ns np nk],is,ip,ik);

    % ---- biophysical constants matching tc_twocmp_stp_all8 --------------
    VL = -70; %#ok<NASGU> % not used by this observer, retained for clarity
    VE =  60; VI = -90; VN = 10; VB = -100;
    mg_switch  = @(V) 1./(1 + 0.28*exp(-0.062*V));
    dmg_switch = @(V) (0.28*0.062*exp(-0.062*V)) ./ (1 + 0.28*exp(-0.062*V)).^2;

    % ---- compartment-placement weights ----------------------------------
    wE_dend_default = [
        0.70 0.80;  % 1 SS
        0.85 0.90;  % 2 SP
        0.45 0.55;  % 3 SI
        0.85 0.90;  % 4 DP
        0.45 0.55;  % 5 DI
        0.80 0.85;  % 6 TP
        0.55 0.65;  % 7 RT
        0.60 0.70]; % 8 RL
    wI_soma_default = 0.90*ones(8,2);

    if isfield(weights,'w_dend') && ~isempty(weights.w_dend)
        wE_dend = weights.w_dend;
    elseif isfield(P,'w_dend') && ~isempty(P.w_dend)
        wE_dend = P.w_dend;
    else
        wE_dend = wE_dend_default;
    end

    if isfield(weights,'w_soma') && ~isempty(weights.w_soma)
        wI_soma = weights.w_soma;
    elseif isfield(P,'w_soma') && ~isempty(P.w_soma)
        wI_soma = P.w_soma;
    else
        wI_soma = wI_soma_default;
    end

    wE_dend = min(max(wE_dend,0),1);
    wI_soma = min(max(wI_soma,0),1);

    % ---- sparse triplets --------------------------------------------------
    rows = []; cols = []; vals = [];

    % optional diagnostic current values at expansion point
    IE = zeros(ns,np); II = zeros(ns,np); Dpop = zeros(ns,np);

    for is = 1:ns
        for p = 1:8
            sw = weights.src(p);
            if sw == 0, continue; end

            Vs = x0(is,p,iVs); Vd = x0(is,p,iVd);
            GE = x0(is,p,iGE); GN = x0(is,p,iGN);
            GI = x0(is,p,iGI); GB = x0(is,p,iGB);

            % compartmentalised conductances at the expansion point
            GE_d = wE_dend(p,1)*GE;
            GN_d = wE_dend(p,2)*GN;
            GI_s = wI_soma(p,1)*GI;
            GB_s = wI_soma(p,2)*GB;

            mgD  = mg_switch(Vd);
            dmgD = dmg_switch(Vd);

            % signed currents, consistent with the membrane equation:
            % excitation tends positive, inhibition tends negative near rest.
            IAMPA = GE_d * (VE - Vd);
            INMDA = GN_d * (VN - Vd) * mgD;
            IGABA = GI_s * (VI - Vs);
            IGB   = GB_s * (VB - Vs);

            IE(is,p) = weights.useAMPA*IAMPA + weights.useNMDA*INMDA;
            II(is,p) = weights.useGABAA*IGABA + weights.useGABAB*IGB;
            Dpop(is,p) = sw * (IE(is,p) - weights.alpha*II(is,p));

            % ---- derivatives wrt states ---------------------------------
            % Excitatory dendritic current derivatives
            if weights.useAMPA
                % d[GE*w*(VE-Vd)] / dGE and dVd
                cGE = sw * wE_dend(p,1) * (VE - Vd);
                cVd = sw * (-GE_d);
                rows = [rows, is, is];
                cols = [cols, idx(is,p,iGE), idx(is,p,iVd)];
                vals = [vals, cGE, cVd];
            end

            if weights.useNMDA
                % d[GN*w*(VN-Vd)*mg(Vd)] / dGN and dVd
                h    = (VN - Vd) * mgD;
                dhdV = -mgD + (VN - Vd) * dmgD;
                cGN  = sw * wE_dend(p,2) * h;
                cVd  = sw * GN_d * dhdV;
                rows = [rows, is, is];
                cols = [cols, idx(is,p,iGN), idx(is,p,iVd)];
                vals = [vals, cGN, cVd];
            end

            % Inhibitory soma current derivatives for -alpha*I_I.
            if weights.useGABAA
                % -alpha * [GI*w*(VI-Vs)]
                cGI = -weights.alpha * sw * wI_soma(p,1) * (VI - Vs);
                cVs =  weights.alpha * sw * GI_s;
                rows = [rows, is, is];
                cols = [cols, idx(is,p,iGI), idx(is,p,iVs)];
                vals = [vals, cGI, cVs];
            end

            if weights.useGABAB
                % -alpha * [GB*w*(VB-Vs)]
                cGB = -weights.alpha * sw * wI_soma(p,2) * (VB - Vs);
                cVs =  weights.alpha * sw * GB_s;
                rows = [rows, is, is];
                cols = [cols, idx(is,p,iGB), idx(is,p,iVs)];
                vals = [vals, cGB, cVs];
            end
        end
    end

    Wsrc = sparse(rows,cols,vals,ns,N);

    if weights.normaliseRows
        rn = sqrt(sum(Wsrc.^2,2));
        rn(rn == 0) = 1;
        Wsrc = spdiags(1./rn,0,ns,ns)*Wsrc;
    end

    if weights.projectLeadfield && isfield(M,'L') && ~isempty(M.L)
        J = M.L * Wsrc;
    else
        J = Wsrc;
    end

    obs = struct();
    obs.kind       = 'RWS_synaptic_current_proxy_linearised';
    obs.alpha      = weights.alpha;
    obs.tauE       = weights.tauE;
    obs.tauI       = weights.tauI;
    obs.src        = weights.src;
    obs.wE_dend    = wE_dend;
    obs.wI_soma    = wI_soma;
    obs.IE         = IE;
    obs.II         = II;
    obs.Dpop       = Dpop;
    obs.note       = ['Static J is the local derivative of I_E,dend - alpha*I_I,soma. ' ...
                      'Apply exp(-s*tauE) to excitatory terms and exp(-s*tauI) to inhibitory terms in frequency-domain code if true delays are required.'];
end
