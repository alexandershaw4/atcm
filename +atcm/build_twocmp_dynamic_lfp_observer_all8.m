function [C0,obs] = build_twocmp_dynamic_lfp_observer_all8(M,P,weights,x0)
% BUILD_TWOCMP_DYNAMIC_LFP_OBSERVER_ALL8
% -------------------------------------------------------------------------
% Dynamic / frequency-domain LFP observation model for tc_twocmp_stp_all8.
%
% This is the dynamic counterpart of build_twocmp_lfp_observer_J_all8.m.
% Instead of collapsing the observation model to one static C/J matrix, it
% returns separate linearised excitatory and inhibitory current observers:
%
%       C_E = d I_E,dend / d x
%       C_I = d I_I,soma / d x
%
% and a frequency/Laplace-domain wrapper:
%
%       C(s) = L * [ exp(-s*tauE) C_E - alpha exp(-s*tauI) C_I ]
%
% where L is M.L when present and requested.  This implements the
% Mazzoni/Linden/Einevoll-style reference weighted synaptic-current proxy
% with a true AMPA/NMDA delay in transfer-function code.
%
% Usage:
%       [C0,obs] = build_twocmp_dynamic_lfp_observer_all8(M,P,weights);
%       Cw       = twocmp_dynamic_lfp_observer_eval(obs,1i*2*pi*f);
%
% The static zero-frequency observer C0 is returned for compatibility:
%       C0 = C(s=0) = L * (C_E - alpha*C_I)
%
% State layout assumed from tc_twocmp_stp_all8:
%   1 Vs, 2 gE, 3 gI, 4 gN, 5 gB, 6 gM, 7 gH, 8 Vd, 9 R, 10 uSTP
%
% Populations:
%   1 SS, 2 SP, 3 SI, 4 DP, 5 DI, 6 TP, 7 RT, 8 RL/RC
% -------------------------------------------------------------------------

    if nargin < 3 || isempty(weights), weights = struct(); end
    if nargin < 4 || isempty(x0), x0 = M.x; end

    % Defaults.  These are deliberately conservative and can be estimated
    % later as observation parameters if desired.
    if ~isfield(weights,'src') || isempty(weights.src)
        weights.src = zeros(1,8);
        weights.src(2) = 1.0;   % superficial pyramidal / main cortical dipole
        weights.src(4) = 0.6;   % deep pyramidal / secondary cortical dipole
    end
    if ~isfield(weights,'alpha') || isempty(weights.alpha), weights.alpha = 1.65; end
    if ~isfield(weights,'tauE')  || isempty(weights.tauE),  weights.tauE  = 0.006; end
    if ~isfield(weights,'tauI')  || isempty(weights.tauI),  weights.tauI  = 0.000; end
    if ~isfield(weights,'useAMPA')  || isempty(weights.useAMPA),  weights.useAMPA  = true; end
    if ~isfield(weights,'useNMDA')  || isempty(weights.useNMDA),  weights.useNMDA  = true; end
    if ~isfield(weights,'useGABAA') || isempty(weights.useGABAA), weights.useGABAA = true; end
    if ~isfield(weights,'useGABAB') || isempty(weights.useGABAB), weights.useGABAB = true; end
    if ~isfield(weights,'projectLeadfield') || isempty(weights.projectLeadfield), weights.projectLeadfield = true; end
    if ~isfield(weights,'normaliseRows') || isempty(weights.normaliseRows), weights.normaliseRows = false; end

    ns = size(M.x,1); np = size(M.x,2); nk = size(M.x,3);
    if np ~= 8 || nk < 8
        error('Expected M.x to have size [ns x 8 x nk] with nk >= 8.');
    end
    x0 = reshape(x0,ns,np,nk);
    N  = numel(M.x);

    iVs = 1; iGE = 2; iGI = 3; iGN = 4; iGB = 5; iVd = 8;
    idx = @(is,ip,ik) sub2ind([ns np nk],is,ip,ik);

    % Reversal potentials and NMDA Mg block, matched to tc_twocmp_stp_all8.
    VE = 60; VI = -90; VN = 10; VB = -100;
    mg_switch  = @(V) 1./(1 + 0.28*exp(-0.062*V));
    dmg_switch = @(V) (0.28*0.062*exp(-0.062*V)) ./ (1 + 0.28*exp(-0.062*V)).^2;

    % Placement defaults: excitatory currents mainly dendritic, inhibition
    % mainly somatic.  These can be overridden by P.w_dend/P.w_soma or weights.
    wE_dend_default = [
        0.70 0.80;  % SS
        0.85 0.90;  % SP
        0.45 0.55;  % SI
        0.85 0.90;  % DP
        0.45 0.55;  % DI
        0.80 0.85;  % TP
        0.55 0.65;  % RT
        0.60 0.70]; % RL
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

    rowsE=[]; colsE=[]; valsE=[];
    rowsI=[]; colsI=[]; valsI=[];
    IE = zeros(ns,np); II = zeros(ns,np);

    for is = 1:ns
        for p = 1:8
            sw = weights.src(p);
            if sw == 0, continue; end

            Vs = x0(is,p,iVs); Vd = x0(is,p,iVd);
            GE = x0(is,p,iGE); GI = x0(is,p,iGI);
            GN = x0(is,p,iGN); GB = x0(is,p,iGB);

            GE_d = wE_dend(p,1)*GE;
            GN_d = wE_dend(p,2)*GN;
            GI_s = wI_soma(p,1)*GI;
            GB_s = wI_soma(p,2)*GB;

            mgD  = mg_switch(Vd);
            dmgD = dmg_switch(Vd);

            % Diagnostic currents at expansion point.
            IAMPA = GE_d * (VE - Vd);
            INMDA = GN_d * (VN - Vd) * mgD;
            IGABA = GI_s * (VI - Vs);
            IGB   = GB_s * (VB - Vs);
            IE(is,p) = weights.useAMPA*IAMPA + weights.useNMDA*INMDA;
            II(is,p) = weights.useGABAA*IGABA + weights.useGABAB*IGB;

            % Excitatory observer derivatives: sw*d(I_AMPA+I_NMDA)/dx.
            if weights.useAMPA
                rowsE = [rowsE, is, is]; %#ok<AGROW>
                colsE = [colsE, idx(is,p,iGE), idx(is,p,iVd)]; %#ok<AGROW>
                valsE = [valsE, sw*wE_dend(p,1)*(VE - Vd), sw*(-GE_d)]; %#ok<AGROW>
            end

            if weights.useNMDA
                h    = (VN - Vd)*mgD;
                dhdV = -mgD + (VN - Vd)*dmgD;
                rowsE = [rowsE, is, is]; %#ok<AGROW>
                colsE = [colsE, idx(is,p,iGN), idx(is,p,iVd)]; %#ok<AGROW>
                valsE = [valsE, sw*wE_dend(p,2)*h, sw*GN_d*dhdV]; %#ok<AGROW>
            end

            % Inhibitory observer derivatives: sw*d(I_GABAA+I_GABAB)/dx.
            % These currents are signed; near rest they are negative.  The
            % dynamic wrapper later forms C_E - alpha*C_I.
            if weights.useGABAA
                rowsI = [rowsI, is, is]; %#ok<AGROW>
                colsI = [colsI, idx(is,p,iGI), idx(is,p,iVs)]; %#ok<AGROW>
                valsI = [valsI, sw*wI_soma(p,1)*(VI - Vs), sw*(-GI_s)]; %#ok<AGROW>
            end

            if weights.useGABAB
                rowsI = [rowsI, is, is]; %#ok<AGROW>
                colsI = [colsI, idx(is,p,iGB), idx(is,p,iVs)]; %#ok<AGROW>
                valsI = [valsI, sw*wI_soma(p,2)*(VB - Vs), sw*(-GB_s)]; %#ok<AGROW>
            end
        end
    end

    CEsrc = sparse(rowsE,colsE,valsE,ns,N);
    CIsrc = sparse(rowsI,colsI,valsI,ns,N);

    if weights.normaliseRows
        Ctmp = CEsrc - weights.alpha*CIsrc;
        rn = sqrt(sum(Ctmp.^2,2));
        rn(rn == 0) = 1;
        S = spdiags(1./rn,0,ns,ns);
        CEsrc = S*CEsrc;
        CIsrc = S*CIsrc;
    end

    if weights.projectLeadfield && isfield(M,'L') && ~isempty(M.L)
        L = M.L;
    else
        L = speye(ns);
    end

    obs = struct();
    obs.kind       = 'dynamic_RWS_synaptic_current_proxy_linearised';
    obs.CEsrc      = CEsrc;
    obs.CIsrc      = CIsrc;
    obs.L          = L;
    obs.alpha      = weights.alpha;
    obs.tauE       = weights.tauE;
    obs.tauI       = weights.tauI;
    obs.src        = weights.src;
    obs.wE_dend    = wE_dend;
    obs.wI_soma    = wI_soma;
    obs.IE         = IE;
    obs.II         = II;
    obs.note       = ['Evaluate with twocmp_dynamic_lfp_observer_eval(obs,s). ' ...
                      'For spectra, use s = 1i*2*pi*f. For a static fallback, use C0.'];

    C0 = atcm.twocmp_dynamic_lfp_observer_eval(obs,0);
end
