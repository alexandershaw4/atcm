function [y,details] = twocmp_lfp_proxy_all8(x,P,M,weights)
% TWOCMP_LFP_PROXY_ALL8
% Direct nonlinear source-level current proxy for tc_twocmp_stp_all8 states.
% This is useful for simulations/time-domain checks. For DCM/transfer-function
% fitting, use build_twocmp_lfp_observer_J_all8.m to obtain a linear observer.

    if nargin < 4 || isempty(weights), weights = struct(); end
    if ~isfield(weights,'src') || isempty(weights.src)
        weights.src = zeros(1,8); weights.src(2)=1; weights.src(4)=0.6;
    end
    if ~isfield(weights,'alpha') || isempty(weights.alpha), weights.alpha = 1.65; end
    if ~isfield(weights,'normalise') || isempty(weights.normalise), weights.normalise = false; end
    if ~isfield(weights,'projectLeadfield') || isempty(weights.projectLeadfield), weights.projectLeadfield = true; end

    ns = size(M.x,1); np = size(M.x,2); nk = size(M.x,3);
    x  = reshape(x,ns,np,nk);
    iVs=1; iGE=2; iGI=3; iGN=4; iGB=5; iVd=8;

    VE=60; VI=-90; VN=10; VB=-100;
    mg_switch = @(V) 1./(1 + 0.28*exp(-0.062*V));

    wE_dend_default = [0.70 0.80;0.85 0.90;0.45 0.55;0.85 0.90;0.45 0.55;0.80 0.85;0.55 0.65;0.60 0.70];
    wI_soma_default = 0.90*ones(8,2);
    if isfield(weights,'w_dend') && ~isempty(weights.w_dend), wE_dend = weights.w_dend;
    elseif isfield(P,'w_dend') && ~isempty(P.w_dend), wE_dend = P.w_dend;
    else, wE_dend = wE_dend_default; end
    if isfield(weights,'w_soma') && ~isempty(weights.w_soma), wI_soma = weights.w_soma;
    elseif isfield(P,'w_soma') && ~isempty(P.w_soma), wI_soma = P.w_soma;
    else, wI_soma = wI_soma_default; end

    ysrc = zeros(ns,1); IE=zeros(ns,np); II=zeros(ns,np);
    for is=1:ns
        for p=1:8
            sw = weights.src(p); if sw==0, continue; end
            Vs=x(is,p,iVs); Vd=x(is,p,iVd);
            GE=x(is,p,iGE); GI=x(is,p,iGI); GN=x(is,p,iGN); GB=x(is,p,iGB);
            IEp = wE_dend(p,1)*GE*(VE - Vd) + wE_dend(p,2)*GN*(VN - Vd)*mg_switch(Vd);
            IIp = wI_soma(p,1)*GI*(VI - Vs) + wI_soma(p,2)*GB*(VB - Vs);
            IE(is,p)=IEp; II(is,p)=IIp;
            ysrc(is) = ysrc(is) + sw*(IEp - weights.alpha*IIp);
        end
    end

    if weights.normalise
        ysrc = ysrc - mean(ysrc);
        s = std(ysrc); if s > 0, ysrc = ysrc./s; end
    end

    if weights.projectLeadfield && isfield(M,'L') && ~isempty(M.L)
        y = M.L * ysrc;
    else
        y = ysrc;
    end

    details = struct('IE',IE,'II',II,'ysrc',ysrc,'alpha',weights.alpha,'src',weights.src);
end
