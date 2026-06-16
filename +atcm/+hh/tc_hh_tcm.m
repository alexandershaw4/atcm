function [f,J,D] = tc_hh_tcm(x,u,P,M)
%TC_HH_TCM Full Hodgkin-Huxley thalamo-cortical mass model.
%
% This is a Hodgkin-Huxley extension of tc_hilge2. It keeps the original
% 8-population thalamo-cortical architecture, intrinsic and extrinsic AMPA,
% NMDA, GABA-A and GABA-B connectivity, thalamo-cortical delays, input routing
% and SPM-compatible calling convention, but replaces the passive membrane
% equation with active HH sodium and potassium currents plus voltage-gated M
% and H currents.
%
% FORMAT
%   [f,J,D] = tc_hh_tcm(x,u,P,M)
%
% Populations, per source
%   1 SS  L4 spiny stellate cells
%   2 SP  L2/3 superficial pyramidal cells
%   3 SI  L2/3 inhibitory interneurons
%   4 DP  L5 deep pyramidal cells
%   5 DI  L5 deep interneurons
%   6 TP  L6 thalamo-cortical projection pyramidal cells
%   7 RT  thalamic reticular cells
%   8 RC  thalamo-cortical relay cells
%
% States, per population
%   1  V      membrane voltage, mV
%   2  gA     AMPA synaptic conductance state
%   3  gG     GABA-A synaptic conductance state
%   4  gN     NMDA synaptic conductance state
%   5  gB     GABA-B synaptic conductance state
%   6  pM     M/Kv7 activation gate, 0..1
%   7  rH     H/HCN activation gate, 0..1
%   8  mNa    fast sodium activation gate, 0..1
%   9  hNa    fast sodium inactivation gate, 0..1
%   10 nK     delayed rectifier potassium activation gate, 0..1
%
% Required M.x size
%   ns x 8 x 10
%
% Compatibility with old P fields
%   Uses P.A, P.AN, P.C, P.H, P.Hn, P.Gb, P.T, P.CV, P.S, P.E, P.D,
%   P.CT, P.TC, P.ID, P.Mh, P.Hh and P.scale_NMDA where present. Missing
%   HH-specific fields use conservative defaults.
%
% Optional HH-specific parameter fields, all log-scaled unless stated
%   P.gNa       1 x 8 or ns x 8, sodium maximal conductance multiplier
%   P.gK        1 x 8 or ns x 8, K delayed-rectifier maximal conductance multiplier
%   P.gL        1 x 8 or ns x 8, leak conductance multiplier
%   P.gM        1 x 8 or ns x 8, M/Kv7 maximal conductance multiplier
%   P.gH        1 x 8 or ns x 8, H/HCN maximal conductance multiplier
%   P.phiHH     scalar, 1 x 8 or ns x 8, HH gate-rate multiplier
%   P.VshiftHH  scalar, 1 x 8 or ns x 8, voltage shift for HH gates, mV, linear
%   P.VhalfM    scalar, 1 x 8 or ns x 8, M half-activation voltage, linear
%   P.kM        scalar, 1 x 8 or ns x 8, M slope, linear
%   P.VhalfH    scalar, 1 x 8 or ns x 8, H half-activation voltage, linear
%   P.kH        scalar, 1 x 8 or ns x 8, H slope, linear
%   P.tauM      scalar, 1 x 8 or ns x 8, M time constant multiplier
%   P.tauH      scalar, 1 x 8 or ns x 8, H time constant multiplier
%   P.NaPop     1 x 8 linear population mask for sodium current
%   P.KPop      1 x 8 linear population mask for delayed K current
%   P.MPop      1 x 8 linear population mask for M current
%   P.HPop      1 x 8 linear population mask for H current
%
% Notes
%   The original TCM already has conductance states for ligand-gated synapses.
%   The HH layer implemented here adds voltage-gated active membrane currents
%   under those synaptic drives. This is not a biophysical single-cell model of
%   every neuron in a population; it is a neural-mass approximation with HH
%   membrane excitability. It will be stiffer than tc_hilge2.
%
% AS2026

% ufun = @(t) 0.5 .* double(t > 0.1 & t < 0.15);
% ufun = @(t) rand
% out = atcm.hh.simulate_hh_tcm_time_domain(DCM.M,DCM.M.pE, ...
%     'method','rk4', ...
%     'T',1.0, ...
%     'dt',5e-5,...
%       'u',ufun);


% -------------------------------------------------------------------------
% Unwrap parameters and dimensions
% -------------------------------------------------------------------------
if isstruct(P) && isfield(P,'p')
    P = P.p;
end

ns = size(M.x,1);
np = size(M.x,2);
nk = size(M.x,3);

if np ~= 8
    error('tc_hh_tcm expects 8 populations, matching tc_hilge2.');
end
if nk < 10
    error('tc_hh_tcm expects M.x to have at least 10 states per population.');
end

x = reshape(x,ns,np,nk);

% Keep gates numerically safe when SPM finite-differences away from the usual range.
V   = x(:,:,1);
gA  = max(x(:,:,2),0);
gG  = max(x(:,:,3),0);
gN  = max(x(:,:,4),0);
gB  = max(x(:,:,5),0);
pM  = clip01(x(:,:,6));
rH  = clip01(x(:,:,7));
mNa = clip01(x(:,:,8));
hNa = clip01(x(:,:,9));
nK  = clip01(x(:,:,10));

% -------------------------------------------------------------------------
% Extrinsic connections and input gains
% -------------------------------------------------------------------------
for k = 1:length(P.A)
    A{k}  = exp(P.A{k}); %#ok<AGROW>
    AN{k} = exp(P.AN{k}); %#ok<AGROW>
end
C = exp(P.C);

for k = 1:numel(A)
    L    = (A{k} > exp(-8)) & (A{k}' > exp(-8));
    A{k} = A{k} ./ (1 + 8*L);
end

% -------------------------------------------------------------------------
% Intrinsic connection strengths and fixed routing masks
% -------------------------------------------------------------------------
G  = exp(full(P.H));
Gn = exp(full(P.Hn));
if isfield(P,'Gb')
    Gb = exp(full(P.Gb));
else
    Gb = G;
end

SA = [1 0 0 0 0;  0 1 0 0 0;  0 1 0 0 0;  0 0 0 0 0; ...
      0 0 0 0 0;  0 0 0 0 0;  0 0 0 0 1;  0 0 0 1 0]/8;
SA(:,[3 4 5]) = 0;

SNMDA = [1 0 0 0 0;  0 1 0 0 0;  0 1 0 0 0;  0 0 0 0 0; ...
         0 0 0 0 0;  0 0 0 0 0;  0 0 0 0 1;  0 0 0 1 0]/8;
SNMDA(:,[3 4 5]) = 0;

GEa = [0 0 0 0 0 2 0 2;
       2 2 0 0 0 0 0 0;
       0 2 0 0 0 0 0 0;
       0 2 0 0 0 0 0 0;
       0 0 0 2 0 0 0 0;
       0 0 0 2 0 0 0 0;
       0 0 0 0 0 0 0 2;
       2 0 0 0 0 2 0 0];

GEn = [0 0 0 0 0 2 0 2;
       2 2 2 0 0 0 0 0;
       0 2 2 0 0 0 0 0;
       0 2 0 0 0 0 0 0;
       0 0 0 2 0 0 0 0;
       0 0 0 2 0 0 0 0;
       0 0 0 0 0 0 0 2;
       2 0 0 0 0 2 0 0];

GIa = [8 0 10 0 0 0 0 0;
       0 18 10 0 0 0 0 0;
       0 0 10 0 0 0 0 0;
       0 0 0 8 6 0 0 0;
       0 0 0 0 14 0 0 0;
       0 0 0 0 6 8 0 0;
       0 0 0 0 0 0 8 0;
       0 0 0 0 0 0 8 8];
GIb = GIa;

if isfield(P,'global')
    GEa = GEa * exp(P.global(1));
    GEn = GEn * exp(P.global(2));
    GIa = GIa * exp(P.global(3));
    GIb = GIb * exp(P.global(4));
end

% -------------------------------------------------------------------------
% Synaptic decay rates, as in tc_hilge2
% -------------------------------------------------------------------------
KE = exp(-P.T(:,1))*1000/2.2;
KI = exp(-P.T(:,2))*1000/5;
KN = exp(-P.T(:,3))*1000/100;
KB = exp(-P.T(:,4))*1000/300;

if isfield(P,'T1')
    KE = KE + P.T1(1);
    KN = KN + P.T1(2);
end

% -------------------------------------------------------------------------
% Reversal potentials and membrane constants, mV and scaled conductance units
% -------------------------------------------------------------------------
VL = -70;      % leak
VE =  60;      % AMPA-like excitation / Na-like excitatory synapse
VI = -90;      % GABA-A / chloride-like inhibition
VN =  10;      % NMDA-like excitation
VB = -100;     % GABA-B-like inhibition
ENa =  55;     % fast sodium
EK  = -90;     % potassium and M/Kv7
EH  = -30;     % H/HCN mixed cation current
VR0 = -52;     % firing transform centre

CV = exp(P.CV) .* [128*3 128 64 128 64 128 64 128] / 1000;
CV = rowmatch(CV,ns,np);

% HH maximal conductances are deliberately conservative relative to classical
% squid axon values because this is a neural mass, not a patch of squid axon.
gNaBar = field_or_default(P,'gNa', log([8 10 4 10 4 8 5 8]), ns, np, true);
gKBar  = field_or_default(P,'gK',  log([3 4 2 4 2 4 3 4]), ns, np, true);
gLBar  = field_or_default(P,'gL',  log(ones(1,8)), ns, np, true);
gMBar  = field_or_default(P,'gM',  log([0.4 0.6 0.3 0.8 0.3 1.4 0.3 1.4]), ns, np, true);
gHBar  = field_or_default(P,'gH',  log([0 0 0 0 0 1.2 0 1.6] + 1e-6), ns, np, true);

NaPop = linear_mask(P,'NaPop',[1 1 0.45 1 0.45 0.8 0.45 0.8],ns,np);
KPop  = linear_mask(P,'KPop', [1 1 0.8  1 0.8  1   0.8  1],ns,np);
MPop  = linear_mask(P,'MPop', [0.4 0.6 0.2 0.8 0.2 1.0 0.1 1.0],ns,np);
HPop  = linear_mask(P,'HPop', [0 0 0 0 0 1 0 1],ns,np);

gNaBar = gNaBar .* NaPop;
gKBar  = gKBar  .* KPop;
gMBar  = gMBar  .* MPop;
gHBar  = gHBar  .* HPop;

phiHH = field_or_default(P,'phiHH',log(0.15),ns,np,true);
Vsh   = field_or_default(P,'VshiftHH',0,ns,np,false);

VhalfM = field_or_default(P,'VhalfM',-35,ns,np,false);
kM     = field_or_default(P,'kM',    10,ns,np,false);
VhalfH = field_or_default(P,'VhalfH',-75,ns,np,false);
kH     = field_or_default(P,'kH',    -7,ns,np,false);
tauM   = field_or_default(P,'tauM',log(160),ns,np,true); % ms
tauH   = field_or_default(P,'tauH',log(100),ns,np,true); % ms

scaleNMDA = 1;
if isfield(P,'scale_NMDA')
    scaleNMDA = exp(P.scale_NMDA);
end

% -------------------------------------------------------------------------
% Firing transform for mass-level synaptic drive
% -------------------------------------------------------------------------
VR = VR0 + P.S;       % or VR0 + exp(P.S), if you want strictly positive threshold shifts
R  = 0.25;            % gentler slope than 2/3 for HH voltages
mfire = 1 ./ (1 + exp(-R .* (V - VR)));

% VR = VR0 + P.S;
% R  = 0.25;
% 
% mV  = 1 ./ (1 + exp(-R .* (V - VR)));
% mNa = x(:,:,8).^3 .* x(:,:,9);     % Na activation × availability
% 
% mfire = mV .* mNa;
% mfire = min(max(mfire,0),1);

% Extrinsic effects, same semantic routing as tc_hilge2
% -------------------------------------------------------------------------
a      = zeros(ns,5);
an     = zeros(ns,5);
a(:,1) = A{1}  * mfire(:,2);
a(:,2) = A{2}  * mfire(:,4);
a(:,3) = A{3}  * mfire(:,6);
a(:,4) = A{4}  * mfire(:,7);
a(:,5) = A{5}  * mfire(:,8);
an(:,1)= AN{1} * mfire(:,2);
an(:,2)= AN{2} * mfire(:,4);
an(:,3)= AN{3} * mfire(:,6);
an(:,4)= AN{4} * mfire(:,7);
an(:,5)= AN{5} * mfire(:,8);

BE = exp(P.E)*0.8;

% -------------------------------------------------------------------------
% Flow
% -------------------------------------------------------------------------
f = zeros(size(x));

for i = 1:ns
    dU = u(:) * C(i,1);

    E     = ( G(:,:,i)  .* GEa) * mfire(i,:)';
    ENMDA = ( Gn(:,:,i) .* GEn) * mfire(i,:)';
    I     = ( G(:,:,i)  .* GIa) * mfire(i,:)';
    IB    = ( Gb(:,:,i) .* GIb) * mfire(i,:)';

    E     = (E     + BE + SA    * a(i,:)')  * 2;
    ENMDA = (ENMDA + BE + SNMDA * an(i,:)') * 2;

    if isfield(P,'endo')
        E(2) = E(2) + 2*exp(P.endo(1));
    end

    if numel(u) > 1
        E(8) = E(8) + dU(2);
        E(2) = E(2) + dU(1);
    else
        input_cell = [8 7];
        E(input_cell) = E(input_cell) + dU;
    end

    if isfield(P,'thi')
        E(8)     = E(8)     + exp(P.thi);
        ENMDA(8) = ENMDA(8) + exp(P.thi);
    end

    Vi = V(i,:);
    mb = local_mag_block(Vi,scaleNMDA);

    % Active HH currents, positive inward current form g * gate * (Erev - V).
    I_leak = gLBar(i,:)                    .* (VL  - Vi);
    I_Na   = gNaBar(i,:) .* mNa(i,:).^3 .* hNa(i,:) .* (ENa - Vi);
    I_K    = gKBar(i,:)  .* nK(i,:).^4              .* (EK  - Vi);
    I_M    = gMBar(i,:)  .* pM(i,:)                 .* (EK  - Vi);
    I_H    = gHBar(i,:)  .* rH(i,:)                 .* (EH  - Vi);

    I_syn  = gA(i,:) .* (VE - Vi) + ...
             gG(i,:) .* (VI - Vi) + ...
             gB(i,:) .* (VB - Vi) + ...
             gN(i,:) .* (VN - Vi) .* mb;

    f(i,:,1) = (I_leak + I_Na + I_K + I_M + I_H + I_syn) ./ CV(i,:);

    % Ligand-gated synaptic conductance states, preserved from tc_hilge2.
    f(i,:,2) = (E'     - x(i,:,2)) .* (KE(i,:)');
    f(i,:,3) = (I'     - x(i,:,3)) .* (KI(i,:)');
    f(i,:,5) = (IB'    - x(i,:,5)) .* (KB(i,:)');
    f(i,:,4) = (ENMDA' - x(i,:,4)) .* (KN(i,:)');

    % HH gates. Rates are per ms in the classical equations, converted to s^-1.
    [am,bm,ah,bh,anr,bnr] = hh_rates(Vi + Vsh(i,:));
    f(i,:,8)  = 1000 * phiHH(i,:) .* (am  .* (1 - mNa(i,:)) - bm  .* mNa(i,:));
    f(i,:,9)  = 1000 * phiHH(i,:) .* (ah  .* (1 - hNa(i,:)) - bh  .* hNa(i,:));
    f(i,:,10) = 1000 * phiHH(i,:) .* (anr .* (1 - nK(i,:))  - bnr .* nK(i,:));

    % M and H gates in steady-state/tau form, converted from ms to s^-1.
    pMinf = 1 ./ (1 + exp(-(Vi - VhalfM(i,:)) ./ max(abs(kM(i,:)),1e-6)));
    rHinf = 1 ./ (1 + exp( (Vi - VhalfH(i,:)) ./ max(abs(kH(i,:)),1e-6)));
    f(i,:,6) = 1000 .* (pMinf - pM(i,:)) ./ tauM(i,:);
    f(i,:,7) = 1000 .* (rHinf - rH(i,:)) ./ tauH(i,:);
end

f = spm_vec(f);

[J,D] = deal([]);

if (nargout < 2 || nargout == 50) && nargin < 5
    return
end

J = spm_cat(spm_diff(M.f,x,u,P,M,1));

if nargout < 3 && nargin < 5
    return
end

% -------------------------------------------------------------------------
% Delays, same scaffolding as tc_hilge2 but automatically expanded to nk=10
% -------------------------------------------------------------------------
CT = 8;
TC = 3;

Tc              = zeros(np,np);
Tc([7 8],1:6)   = CT * exp(P.CT);
Tc(1:6,[7 8])   = TC * exp(P.TC);
Tc              = Tc / 1000;
Tc              = kron(ones(nk,nk),kron(Tc,ones(ns,ns)));

ID = [2 1 1 1 1 2 1 2];
ID = ID .* exp(P.ID) / 1000;
ID = repmat(ID,[1 nk]);
ID = repmat(ID(:)',[np*nk,1]);
ID = kron(ID,ones(ns,ns));

D = Tc + ID;

end

% =========================================================================
% Local utilities
% =========================================================================
function y = clip01(x)
y = min(max(x,0),1);
end

function X = rowmatch(x,ns,np)
% Convert scalar, 1 x np, ns x np or vector to ns x np.
if isscalar(x)
    X = repmat(x,ns,np);
elseif isequal(size(x),[1 np])
    X = repmat(x,ns,1);
elseif isequal(size(x),[ns np])
    X = x;
elseif numel(x) == np
    X = repmat(reshape(x,1,np),ns,1);
else
    error('Cannot expand parameter with size %s to ns x np.',mat2str(size(x)));
end
end

function X = field_or_default(P,field,default,ns,np,islog)
if isfield(P,field)
    val = P.(field);
else
    val = default;
end
X = rowmatch(val,ns,np);
if islog
    X = exp(X);
end
end

function Msk = linear_mask(P,field,default,ns,np)
if isfield(P,field)
    Msk = rowmatch(P.(field),ns,np);
else
    Msk = rowmatch(default,ns,np);
end
end

function mb = local_mag_block(vrow,scaleNMDA)
mb = 1 ./ (1 + 0.2 .* exp(-0.062 .* scaleNMDA .* vrow));
end

function [am,bm,ah,bh,an,bn] = hh_rates(V)
% Classical HH rates in ms^-1, using the common modern V in mV convention
% where resting voltage is approximately -65 mV.
am = 0.1  .* vtrap(V + 40,10);
bm = 4.0  .* exp(-(V + 65)./18);
ah = 0.07 .* exp(-(V + 65)./20);
bh = 1.0  ./ (1 + exp(-(V + 35)./10));
an = 0.01 .* vtrap(V + 55,10);
bn = 0.125.* exp(-(V + 65)./80);
end

function y = vtrap(x,yden)
% Numerically stable x/(1-exp(-x/yden)).
z = x ./ yden;
y = x ./ (1 - exp(-z));
small = abs(z) < 1e-6;
y(small) = yden .* (1 + z(small)./2);
end
