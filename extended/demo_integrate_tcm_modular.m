function demo_integrate_tcm_modular()
% DEMO_INTEGRATE_TCM_MODULAR  Numerical simulation of the extended TCM (time-domain).
%
% This uses your compiled mega-model (tcm_modular) with a fixed-step RK4 integrator.

% --- Registries & compile ---
ns  = 1;
Reg = build_default_registries();
Comp = compile_tcm(Reg, ns);

% --- Model struct M ---
np = numel(Reg.Pops);
nk = 1 + numel(Reg.Channels);
M = struct();
M.x = zeros(ns, np, nk);       % state template
M.Comp = Comp;
M.Reg  = Reg;
M.f    = @tcm_modular;         % your state equation: [f,J,D] = tcm_modular(x,u,P,M)

% --- Priors / parameters ---
opts = struct();
opts.free_mode = 'class-tied'; % or 'per-receptor', etc.
Pri  = setup_priors_extended(opts);

M.labels = Pri.labels;
P  = Pri.vec_to_struct(Pri.pE);
PC = Pri.vec_to_struct(Pri.pC); %#ok<NASGU> % (kept for compatibility if you need it)
% Optional: small jitter to break symmetry
%P = jitter_some_params(P, M.labels, 0.05);

% --- Initial condition (vectorised) ---
x0 = spm_vec(M.x);

% --- Input function u(t) ---
% Example: two inputs (as your previous [0;0]) — here we inject a short pulse onto RC
uRC = @(t) (t>1 & t<1.2) * 1.0;  % 200 ms pulse at 1 s
ufun = @(t) [0; uRC(t)];         % adapt length to however many inputs your model expects

% --- Observation (mean of voltages by default) ---
obs = struct(); % leave empty to auto-define mean(V); or set obs.C manually

% --- Integrate ---
dt = 1/600; 1e-3;   % 1 ms
T  = 1.0;    % 5 seconds
[t, X, Y, Meta] = integrate_tcm_modular(P, M, dt, T, ufun, x0, obs);

% --- Unpack and plot voltages ---
V = extract_voltages(X, Meta);  % [N x np], voltage trace per population
figure('Color','w');
plot(t, V, 'LineWidth',1.0);
grid on
xlabel('Time (s)');
ylabel('Membrane potential (a.u.)');
legend(Reg.Pops, 'Location','eastoutside');
title('Extended TCM: membrane potentials');

% --- Also plot the observation (if computed) ---
if ~isempty(Y)
    figure('Color','w');
    plot(t, Y, 'LineWidth',1.2);
    grid on
    xlabel('Time (s)');
    ylabel('Observed signal (a.u.)');
    title('Observation y(t)');
end

end

% ===== helper: extract voltage traces (state #1 for each population) =====
function V = extract_voltages(X, Meta)
N  = size(X,1);
np = Meta.np; nk = Meta.nk;
V  = zeros(N, np);
for p = 1:np
    idx = Meta.index(1, p, 1);      % (s=1, pop=p, state k=1 => V)
    V(:, p) = X(:, idx);
end
end

% ===== optional: tiny jitter on selected parameters =====
function Pout = jitter_some_params(Pin, labels, sd)
Pout = Pin;
% Example: jitter all synaptic weights (labels starting 'w_')
mask = startsWith(labels, 'w_');
for k = find(mask(:).')
    f = matlab.lang.makeValidName(labels{k});
    if isfield(Pout, f)
        Pout.(f) = Pout.(f) + sd*randn;
    end
end
end
