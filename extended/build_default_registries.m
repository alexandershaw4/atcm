% TCM Modular Refactor — Starter Pack (MATLAB)
% This canvas contains a minimal, declarative framework to extend your TCM with
% additional receptors/modulators and non–hard-coded state–state interactions.
% 
% Files included below 
%   1) build_default_registries.m     — Channels/Mods/Pops/Syn/Delays registries
%   2) compile_tcm.m                  — Compile registries into sparse operators
%   3) tcm_modular.m                  — Main state-equation using compiled ops
%   4) modulation_utils.m             — Presyn/postsyn modulation helpers
%   5) stp_utils.m                    — Tsodyks–Markram short-term plasticity
%   6) gates_and_transforms.m         — Nonlinear gates (e.g., NMDA Mg block)
%   7) interaction_kernels.m          — Optional state–state mixing kernel K
%   8) demo_run_tcm_modular.m         — Tiny demo that mirrors your current wiring
%
% You can split these into separate .m files one-to-one. They’re written to be
% drop-in alongside your SPM/DCM toolchain.

% | Population                         | AMPA | NMDA | KA | GABAA | GABAB | M | HCN | SK | NaP | 5HT2A (post) | 5HT1A (post) | ACh M1 (post) | DA D1 (post) | NE β (post) | CB1 (pre) | MOR (pre) |
% | ---------------------------------- | ---- | ---- | -- | ----- | ----- | - | --- | -- | --- | ------------ | ------------ | ------------- | ------------ | ----------- | --------- | --------- |
% | **SS** (L4 spiny stellate)         | ✓    | ✓    | ✓  | ✓     |       |   |     |    |     |              |              |               |              |             |           |           |
% | **SP** (L2/3 superficial pyramids) | ✓    | ✓    | ✓  | ✓     | ✓     | ✓ |     | ✓  | ✓   | ✓            | ✓            | ✓             | ✓            | ✓           | ✓         |           |
% | **SI** (L2/3 interneurons)         | ✓    | ✓    |    | ✓     | ✓     |   |     |    |     |              |              | ✓             |              |             |           | ✓         |
% | **DP** (L5 deep pyramids)          | ✓    | ✓    | ✓  | ✓     | ✓     | ✓ |     | ✓  | ✓   | ✓            | ✓            | ✓             | ✓            | ✓           |           |           |
% | **DI** (L5 deep interneurons)      | ✓    | ✓    |    | ✓     | ✓     |   |     |    |     |              |              | ✓             |              |             |           | ✓         |
% | **TP** (L6 thal. projection)       | ✓    | ✓    | ✓  | ✓     | ✓     | ✓ | ✓   | ✓  | ✓   |              |              |               |              | ✓           |           |           |
% | **RT** (thalamic reticular)        | ✓    | ✓    |    | ✓     | ✓     |   |     |    |     |              |              |               |              |             |           |           |
% | **RC** (thalamic relay)            | ✓    | ✓    | ✓  | ✓     | ✓     | ✓ | ✓   | ✓  | ✓   |              |              |               |              | ✓           |           |           |



% 1) build_default_registries.m
function Reg = build_default_registries()
%BUILD_DEFAULT_REGISTRIES  Declarative registries for a typed-graph TCM.
% Returns a struct Reg with fields:
%  Reg.Pops, Reg.Channels, Reg.Mods, Reg.Syn, Reg.Delays, Reg.STP, Reg.Params

% --- Populations (id order must match your DCM/M.x layout) ---
Reg.Pops = struct('name',{ 'SS','SP','SI','DP','DI','TP','RT','RC' });
np = numel(Reg.Pops);

% --- Channels (conductance states) ---
% name, reversal (mV), tau_name (maps to parameter), gate_id (0=none,1=NMDA Mg), sign (+ depol / - hyper)
% Reg.Channels = struct( ...
%   'name',    { 'AMPA','NMDA','GABAA','GABAB','M','HCN','KA','SK','NaP' }, ...
%   'Erev',    [   60,    10,    -90,   -100,  -52,  -30,   60,  -90,   55 ], ...
%   'tau_name',{'KE','KN','KI','KB','KM','KH','KKA','KSK','KNaP'}, ...
%   'gate_id', [    0,     1,      0,      0,    0,    0,    0,    0,    0 ], ...
%   'sign',    [   +1,    +1,     -1,     -1,   -1,   +1,   +1,   -1,   +1 ]);


Reg.Channels = [
  struct('name','AMPA',  'Erev',  60,  'tau_name','KE',   'gate_id',0, 'sign',+1)
  struct('name','NMDA',  'Erev',  10,  'tau_name','KN',   'gate_id',1, 'sign',+1)
  struct('name','GABAA', 'Erev', -90,  'tau_name','KI',   'gate_id',0, 'sign',-1)
  struct('name','GABAB', 'Erev',-100,  'tau_name','KB',   'gate_id',0, 'sign',-1)
  struct('name','M',     'Erev', -52,  'tau_name','KM',   'gate_id',0, 'sign',-1)
  struct('name','HCN',   'Erev', -30,  'tau_name','KH',   'gate_id',0, 'sign',+1)
  struct('name','KA',    'Erev',  60,  'tau_name','KKA',  'gate_id',0, 'sign',+1)
  struct('name','SK',    'Erev', -90,  'tau_name','KSK',  'gate_id',0, 'sign',-1)
  struct('name','NaP',   'Erev',  55,  'tau_name','KNaP', 'gate_id',0, 'sign',+1)
];

% --- Modulators (slow fields) ---
Reg.Mods = struct( ...
  'name', {'s5HT2A','s5HT1A','ACh_M1','DA_D1','NE_beta','CB1','MOR'}, ...
  'kin_name',{'k_5HT2A','k_5HT1A','k_M1','k_D1','k_beta','k_CB1','k_Mu'} );

% --- Delays (families) ---
Reg.Delays = struct( ...
  'name', {'Fwd','Bwd','Cx2Th','Th2Cx','Intra'}, ...
  'value_ms',{ 8,     8,     60,      20,      2 } );

% --- Synapses (edge list) ---
% Columns: pre post receptor tau_family delay_family stp_id pre_mods post_mods param_name
% Syn = {
% %  pre  post  receptor  delay   stp  pre_mods           post_mods        param
%   'SP', 'SS', 'AMPA',   'Fwd',   0,  {'CB1'},           {},              'w_SP_SS_AMPA';
%   'SP', 'SS', 'NMDA',   'Fwd',   0,  {},                {'s5HT2A'},       'w_SP_SS_NMDA';
%   'DP', 'SP', 'AMPA',   'Bwd',   0,  {},                {},              'w_DP_SP_AMPA';
%   'DP', 'SP', 'NMDA',   'Bwd',   0,  {},                {'s5HT2A','DA_D1'},'w_DP_SP_NMDA';
%   'SP', 'SI', 'AMPA',   'Fwd',   1,  {'CB1','MOR'},     {},              'w_SP_SI_AMPA';
%   'SI', 'DP', 'GABAA',  'Intra', 0,  {},                {},              'w_SI_DP_GABAA';
%   'TP', 'RC', 'AMPA',   'Cx2Th', 0,  {},                {},              'w_TP_RC_AMPA';
%   'RC', 'RT', 'AMPA',   'Th2Cx', 0,  {},                {},              'w_RC_RT_AMPA';
%   'RT', 'RC', 'GABAA',  'Th2Cx', 0,  {},                {},              'w_RT_RC_GABAA';
% };

pop = {'SS','SP','SI','DP','DI','TP','RT','RC'};  % must match your M.x order

M = [ ...
  0.1250 0      0      0      0      0      0      0.1250;  % to SS
  0.1250 0.1250 0.1250 0      0      0      0      0     ;  % to SP
  0      0.1250 0.1250 0      0      0      0      0     ;  % to SI
  0      0.1250 0      0.1250 0.1250 0      0      0     ;  % to DP
  0      0      0      0      0.1250 0      0      0     ;  % to DI
  0      0      0      0.1250 0.1250 0.1250 0      0     ;  % to TP
  0      0      0      0      0      0      0      0     ;  % to RT
  0.1250 0      0      0      0      0.1250 0.1250 0     ]; % to RC

Syn = syn_from_matrix(M, pop);

Reg.Syn = Syn;%cell2table(Syn, 'VariableNames', ...
  %{'pre','post','receptor','delay_family','stp_id','pre_mods','post_mods','param'});

% --- Short-term plasticity defaults (id -> params) ---
Reg.STP = struct('U0',0.2,'D',0.4,'F',0.2); % global defaults; can extend per-edge

% --- Parameter names / priors (skeleton) ---
Reg.Params = struct();
Reg.Params.gain_mod_pre  = 0.0;  % log-sensitivity default
Reg.Params.gain_mod_post = 0.0;
Reg.Params.tau = struct('KE',2.2,'KI',6,'KN',100,'KB',300,'KM',160,'KH',100,'KKA',4,'KSK',6,'KNaP',10); % ms
Reg.Params.capacitance = 128*ones(1,np); % pF scaled; adjust per-pop
Reg.Params.Vrev = struct('VL',-70,'VE',60,'VI',-90,'VN',10,'VB',-100,'VM',-52,'VH',-30);
end

