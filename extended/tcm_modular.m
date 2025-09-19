function [f,J,D] = tcm_modular(x,u,P,M)
%TCM_MODULAR  State equations using declarative compiled operators.
% x reshaped as (ns,np,nk): states: [V, g(AMPA), g(NMDA), g(GABAA), g(GABAB), g(M), g(HCN), ...]
% P carries log-params for edge weights, taus, mod sensitivities, etc.

Comp = M.Comp;
Reg = M.Reg;


% --- Dimensions & reshape ---
ns = size(M.x,1); np = size(M.x,2); nk = size(M.x,3);
x  = reshape(x,ns,np,nk);

% --- Firing rates ---
VR = -52; R = 2/3; V = x(:,:,1);
m  = 1./(1 + exp(-R*(V - VR))); % sigmoid; replace with your spm_Ncdf flavour if preferred

% --- Exogenous input (optional: to selected pops) ---
Ex = zeros(ns,np);
if ~isempty(u)
  Ex(:,2) = Ex(:,2) + u(1); % SP
  if numel(u)>1, Ex(:,8) = Ex(:,8) + u(2); end % RC
end

% --- Build per-receptor drive (rate-based) ---
% Edge weights come from P.(param) in Syn table. For brevity, use exp() lookup.
receptors = fieldnames(Comp.A);
Drive = struct();
for r = 1:numel(receptors)
  rec = receptors{r};
  A = Comp.A.(rec);
  w = edge_weights(Reg.Syn, P);             % vector length = Nedges (uniform for now)
  % Scatter edge weights into A’s nonzeros (simple scalar here); extend to per-edge if storing separately
  % For now treat A as 0/1 topology and fold weights into a diag on pres -> fast multiply via (A*m)
  m_vec = reshape(m', [], ns); m_vec = m_vec(:); % stack per source
  d = A * m_vec;                                 % post-syn drive
  Drive.(rec) = reshape(d, np, ns)';             % back to (ns,np)
end

% --- Add background and exogenous per receptor family ---
if isfield(P,'BE'), BE = exp(P.BE); else, BE = 0.8; end
for r = 1:numel(receptors)
  Drive.(receptors{r}) = Drive.(receptors{r}) + BE + Ex; %#ok<AGROW>
end

% --- Apply presynaptic modulation (CB1, MOR, etc.) ---
Drive = apply_presyn_modulation(Drive, Reg, Comp, P);

% --- Conductance ODEs per channel ---
Erev = containers.Map({Reg.Channels.name}, num2cell([Reg.Channels.Erev]));
Tau  = get_channel_taus(Reg, P, ns, np);
G    = struct();               % conductances (for voltage eq)
f    = zeros(ns, np, nk);      % <-- derivatives container (do NOT set f = x)

% Cache fields some channels need
V     = x(:,:,1);
nmda_idx = find(strcmp({Reg.Channels.name},'NMDA'), 1);
gNMDA = x(:,:, 1 + nmda_idx);

for c = 1:numel(Reg.Channels)
  ch = Reg.Channels(c).name;

  if isfield(Drive, ch)
    target = Drive.(ch);                             % synaptic receptors
  else
    target = intrinsic_target(ch, V, gNMDA, P);      % intrinsic channels (M, HCN, SK, NaP)
  end

  % postsynaptic modulation (e.g., 5HT2A)
  gain_post = apply_postsyn_modulation(ch, Reg, Comp, P, ns, np);
  if isempty(gain_post), gain_post = 1; end

  g  = x(:,:, 1 + c);                % current conductance
  dg = (gain_post .* target - g) .* Tau.(ch);
  f(:,:, 1 + c) = dg;                % <-- write derivative (do NOT update x here)
  G.(ch)        = g;                 % cache for voltage
end

% --- Voltage ODE ---
VL = -70; GL = 1; Ccap = 128/1000;   % replace with per-pop if you have them
sum_curr = GL*(VL - V);
for c = 1:numel(Reg.Channels)
  ch = Reg.Channels(c).name; Er = Erev(ch);
  gate = gate_value(ch, V, P);       % NMDA Mg block etc.
  sum_curr = sum_curr + G.(ch) .* (Er - V) .* gate;
end
% Optional state–state interaction kernel
if isfield(P,'K')
  sum_curr = sum_curr + apply_interaction_kernel(V, x, P, Reg);
end

f(:,:,1) = sum_curr ./ Ccap;         % <-- voltage derivative

% Vectorise & outputs
f = spm_vec(f);
J = [];
D = kron(ones(nk,nk),kron(Comp.delay,ones(ns,ns)));                      
if nargout > 1
  J = spm_cat(spm_diff(M.f, x, u, P, M, 1));
end
end

function w = edge_weights(Syn, P)
% Return vector of edge weights (exp of named params); simple version
w = ones(height(Syn),1);
% Example: if P holds fields with names in Syn.param (log-space)
for k = 1:height(Syn)
  nm = Syn.param{k};
  if isfield(P, nm), w(k) = exp(P.(nm)); end
end
end

function Tau = get_channel_taus(Reg, P, ns, np)
Tau = struct();
for c = 1:numel(Reg.Channels)
  nm = Reg.Channels(c).tau_name;
  val = Reg.Params.tau.(nm); % default (ms)
  if isfield(P, nm), val = exp(P.(nm)); end
  Tau.(Reg.Channels(c).name) = (1./(val/1000)) * ones(ns,np); % rate = 1/tau(s)
end
end
