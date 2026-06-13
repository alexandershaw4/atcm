function feat = extract_ngnmm_features(P, M, x)
% extract_ngnmm_features
% Extract interpretable dynamical features from fitted NGNMM parameters
% and states.
%
% Inputs
%   P : parameter structure, e.g. DCM.Ep or DCM.M.pE
%   M : model structure, e.g. DCM.M
%   x : state array, usually DCM.M.x or fixed point x_fp
%
% Output
%   feat : structure of interpretable NGNMM summaries
%
% AS20206

    if nargin < 3 || isempty(x)
        x = M.x;
    end

    [f,J,D,aux] = feval(M.f,x,0,P,M);

    R = aux.R;
    V = aux.V;
    g = aux.g;
    h = aux.h;

    Kappa = aux.Kappa;
    Erev  = aux.Erev;
    Eeff  = aux.Eeff;

    alpha = exp(P.logAlpha(:))';
    tau_alpha = 1 ./ alpha;

    Delta = exp(P.logDelta(:))';
    eta   = P.eta(:)';

    % Excitatory / inhibitory source split
    exc_sources = M.popType > 0;
    inh_sources = M.popType < 0;

    Gexc = zeros(size(R));
    Ginh = zeros(size(R));

    for s = 1:M.Ns
        for p = 1:M.Npop
            Gexc(s,p) = sum(Kappa(p,exc_sources) .* R(s,exc_sources));
            Ginh(s,p) = sum(Kappa(p,inh_sources) .* R(s,inh_sources));
        end
    end

    EI_conductance_ratio = Gexc ./ max(eps,Ginh);

    % Eigenmode diagnostics
    ev = eig(J);
    Hz = abs(imag(ev)) ./ (2*pi);
    damping = real(ev);
    Q = abs(imag(ev)) ./ max(eps,-2*damping);

    % Observation weights
    if isfield(P,'J')
        Obs = reshape(exp(P.J(:)),size(M.x));
    else
        Obs = [];
    end

    % Package
    feat = struct();

    feat.R = R;
    feat.V = V;
    feat.g = g;
    feat.h = h;

    feat.meanR = squeeze(mean(R,1));
    feat.meanV = squeeze(mean(V,1));
    feat.meanG = squeeze(mean(g,1));
    feat.meanH = squeeze(mean(h,1));

    feat.Kappa = Kappa;
    feat.Erev = Erev;
    feat.Eeff = Eeff;
    feat.drive_g = aux.drive_g;

    feat.Gexc = Gexc;
    feat.Ginh = Ginh;
    feat.EI_conductance_ratio = EI_conductance_ratio;

    feat.alpha = alpha;
    feat.tau_alpha = tau_alpha;
    feat.Delta = Delta;
    feat.eta = eta;

    feat.sync = aux.sync;
    feat.phase = aux.phase;

    feat.Jacobian = J;
    feat.DelayMatrix = D;
    feat.eigenvalues = ev;
    feat.modeHz = Hz;
    feat.modeDamping = damping;
    feat.modeQ = Q;

    feat.Observer = Obs;

    % Some named loop gains for the 8-population model
    SS = 1; SP = 2; SI = 3; DP = 4; DI = 5; TP = 6; RT = 7; RL = 8;

    loops = struct();
    loops.superficial_EI = Kappa(SI,SS) * Kappa(SP,SI);
    loops.deep_EI        = Kappa(DI,DP) * Kappa(DP,DI);
    loops.RT_RL          = Kappa(RL,RT) * Kappa(RT,RL);
    loops.TP_RT          = Kappa(RT,TP) * Kappa(TP,RT);
    loops.RL_cortex      = Kappa(SS,RL) + Kappa(SP,RL) + Kappa(SI,RL);
    loops.cortex_thal    = Kappa(TP,DP) + Kappa(RT,DP) + Kappa(RL,DP);

    feat.loops = loops;

    % Convenient population table
    pop = M.pop(:);

    feat.populationTable = table( ...
        pop, ...
        feat.meanR(:), ...
        feat.meanV(:), ...
        feat.meanG(:), ...
        alpha(:), ...
        tau_alpha(:), ...
        eta(:), ...
        Delta(:), ...
        Eeff(:), ...
        squeeze(mean(Gexc,1))', ...
        squeeze(mean(Ginh,1))', ...
        squeeze(mean(EI_conductance_ratio,1))', ...
        squeeze(mean(aux.sync,1))', ...
        'VariableNames', { ...
            'Population', ...
            'MeanRate_R', ...
            'MeanVoltage_V', ...
            'MeanConductance_g', ...
            'AlphaRate', ...
            'AlphaTau', ...
            'Eta', ...
            'Delta', ...
            'Eeff', ...
            'Gexc', ...
            'Ginh', ...
            'EIConductanceRatio', ...
            'Synchrony'});
end