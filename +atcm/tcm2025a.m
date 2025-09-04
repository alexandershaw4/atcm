function [f,J,D] = tcm2025a(x,u,P,M)
% State equations for an extended canonical thalamo-cortical neural-mass model.
%
% Conductance-based neural mass with AMPA, NMDA, GABA-A, GABA-B (+ optional M & H).
% Populations (per source): 
% 1 SS (L4 spiny stellate), 2 SP (L2/3 sup. pyramids), 3 SI (L2/3 interneurons),
% 4 DP (L5 deep pyramids),   5 DI (L5 deep interneurons),
% 6 TP (L6 thal proj pyramids), 7 RT (reticular), 8 RC (relay)
%
% Outputs:
%   f : vectorised state derivatives (spm_vec order)
%   J : Jacobian (df/dx) if requested
%   D : delay matrix if requested
%
% Dr Alexander Shaw | 2020 | alexandershaw4[@]gmail.com

% Allow P.p struct wrapping
if isstruct(P) && isfield(P,'p')
    P = P.p;
end

% --- Toggle optional M/H channels on L6 (TP) & thalamic relay (RC)
IncludeMH = 1;

% --- Dimensions and reshape state ----------------------------------------
ns = size(M.x,1);              % sources
np = size(M.x,2);              % populations
nk = size(M.x,3);              % states per population
x  = reshape(x,ns,np,nk);      % hidden states

% --- inter-state effects ----------------------------
Xpop   = 0.1 * tanh(full(P.Xpop));    % 8x8 population couplings (small)
Xstate = 0.1 * tanh(full(P.Xstate));  % 7x7 state-type couplings (small)
Xgain  = exp(P.Xgain);                % scalar gain (start ~ e^-3..e^-2)

% Masks to control biology (edit as you like)
Mstate = ones(nk);               % allow all by default...
Mstate(1, :) = 0;               % ...except no contributions INTO V

% Optional: forbid self-population short-cuts or allow only within-source
%Mpop   = ones(np);               % allow cross-pop too\
Mpop = eye(8);                % uncomment to keep within-pop only

Kstate = Xstate .* Mstate;      % 7x7
Kpop   = Xpop   .* Mpop;        % 8x8

% Build block operator once (outside the ns loop)
Kblock = kron(Kpop, Kstate);    % (8*7) x (8*7) = 56x56


% --- Extrinsic connections (A / AN) and input gains -----------------------
for i = 1:length( P.A )
    A{i}  = exp(P.A{i});
    AN{i} = exp(P.AN{i});
end
C  = exp(P.C);

% Damp reciprocal (lateral) strengths on A (as in original)
for i = 1:numel(A)
    L    = (A{i} > exp(-8)) & (A{i}' > exp(-8));
    A{i} = A{i}./(1 + 8*L);
end

% --- Intrinsic (within-source) connection strengths -----------------------
G  = exp(full(P.H));
Gn = exp(full(P.Hn));

if isfield(P,'Gb')
    Gb = exp(full(P.Gb));
else
    Gb = G;
end

% --- Extrinsic routing masks (fixed switches as per original) -------------
% AMPA-mediated
SA = [1 0 0 0 0;  0 1 0 0 0;  0 1 0 0 0;  0 0 0 0 0; ...
      0 0 0 0 0;  0 0 0 0 0;  0 0 0 0 1;  0 0 0 1 0]/8;
SA(:,[3 4 5]) = 0; % ketamine study restriction

% NMDA-mediated
SNMDA = [1 0 0 0 0;  0 1 0 0 0;  0 1 0 0 0;  0 0 0 0 0; ...
         0 0 0 0 0;  0 0 0 0 0;  0 0 0 0 1;  0 0 0 1 0]/8;
SNMDA(:,[3 4 5]) = 0; % ketamine study restriction

% --- Intrinsic switches (fixed topology) ----------------------------------
% Excitatory (AMPA/NMDA)
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

% Inhibitory (GABA-A/B)
GIa = [8 0 10 0 0 0 0 0;
       0 18 10 0 0 0 0 0;
       0 0 10 0 0 0 0 0;
       0 0 0 8 6 0 0 0;
       0 0 0 0 14 0 0 0;
       0 0 0 0 6 8 0 0;
       0 0 0 0 0 0 8 0;
       0 0 0 0 0 0 8 8];
GIb = GIa;

% --- Channel time constants (decay rates) ---------------------------------
KE = exp(-P.T(:,1))*1000/2.2;   % AMPA
KI = exp(-P.T(:,2))*1000/5;     % GABA-A
KN = exp(-P.T(:,3))*1000/100;   % NMDA
KB = exp(-P.T(:,4))*1000/300;   % GABA-B

% Trial-specific (optional) AMPA/NMDA mod for LTP task
if isfield(P,'T1')
    KE = KE + P.T1(1);
    KN = KN + P.T1(2);
end

% --- Reversal potentials & optional M/H channels --------------------------
VL = -70; VE =  60; VI = -90; VR = -52; VN = 10; VB = -100;

if IncludeMH
    VM = -52; VH = -30;               % M / H
    GIm = diag(4*[1 1 1 1 1 1 1 1].*exp(P.Mh(:)'));
    GIh = diag(4*[0 0 0 0 0 1 0 1].*exp(P.Hh(:)'));
    KM  = exp(-P.T(:,5))*1000/160;    % M current rate
    KH  = exp(-P.T(:,6))*1000/100;    % H current rate
end

% --- Membrane capacitances & leak -----------------------------------------
CV = exp(P.CV).*[128*3 128 64 128 64 128 64 128]/1000;  % per population
GL = 1;

% --- Mean-field excitability shifts & sigmoid firing ----------------------
VR = VR + exp(P.S);        % slope shift
R  = 2/3;
FF = 1./(1 + exp(-R.*(x(:,:,1) - VR)));
RS = 30;
FF(x(:,:,1) >= VR) = 1;
FF(x(:,:,1) >= RS) = 0;
m  = FF;                   % firing proxy

% H-channel mean firing surrogate (match original form)
if IncludeMH
    h = 1 - 1./(1 + exp(-(2/3).*(x(:,:,1) - VH)));
end

% --- Extrinsic effects (per source) ---------------------------------------
a      = zeros(ns,5);
an     = zeros(ns,5);
a(:,1) = A{1} * m(:,2);   % SP->SS (F)
a(:,2) = A{2} * m(:,4);   % DP->SP (B)
a(:,3) = A{3} * m(:,6);   % TP
a(:,4) = A{4} * m(:,7);   % RT
a(:,5) = A{5} * m(:,8);   % RC
an(:,1)= AN{1} * m(:,2);
an(:,2)= AN{2} * m(:,4);
an(:,3)= AN{3} * m(:,6);
an(:,4)= AN{4} * m(:,7);
an(:,5)= AN{5} * m(:,8);

% --- Background drive ------------------------------------------------------
BE = exp(P.E)*0.8;

% --- Optional global scaling of intrinsic blocks --------------------------
if isfield(P,'global')
    GEa = GEa * exp(P.global(1));
    GEn = GEn * exp(P.global(2));
    GIa = GIa * exp(P.global(3));
    GIb = GIb * exp(P.global(4));
end

% --- Flow over sources/populations ----------------------------------------
f = x;   % preallocate like x

for i = 1:ns
    % Input scaling (per source)
    dU = u(:) * C(i,1);

    % Intrinsic coupling (parameterised)
    E     = ( G(:,:,i).*GEa ) * m(i,:)';   % AMPA
    ENMDA = ( Gn(:,:,i).*GEn ) * m(i,:)';  % NMDA
    I     = ( G(:,:,i).*GIa ) * m(i,:)';   % GABA-A
    IB    = ( Gb(:,:,i).*GIb ) * m(i,:)';  % GABA-B

    % Optional intrinsic M/H (non-parameterised) currents
    if IncludeMH
        Im = GIm * m(i,:)';
        Ih = GIh * h(i,:)';
    end

    % Extrinsic excitatory + background (×2 as in original)
    E     = (E     + BE + SA   * a (i,:)') * 2;
    ENMDA = (ENMDA + BE + SNMDA* an(i,:)') * 2;

    % Optional endogenous boost to SP
    if isfield(P,'endo')
        E(2) = E(2) + 2*exp(P.endo(1));
    end

    % Exogenous input routing
    if numel(u) > 1
        E(8) = E(8) + dU(2);  % relay
        E(2) = E(2) + dU(1);  % superficial pyramids
    else
        input_cell  = [8 7];  % thalamus (relay & reticular)
        E(input_cell) = E(input_cell) + dU;
    end

    % Direct thalamic current (optional)
    if isfield(P,'thi')
        E(8)     = E(8)     + exp(P.thi);
        ENMDA(8) = ENMDA(8) + exp(P.thi);
    end

    % --- Voltage equations -------------------------------------------------
    if ~IncludeMH
        % Use external mg_switch (as in original)
        f(i,:,1) = ( GL*(VL - x(i,:,1)) + ...
                     1.0*x(i,:,2).*(VE - x(i,:,1)) + ...
                     1.0*x(i,:,3).*(VI - x(i,:,1)) + ...
                     1.0*x(i,:,5).*(VB - x(i,:,1)) + ...
                     1.0*x(i,:,4).*(VN - x(i,:,1)).*mg_switch(x(i,:,1)) ) ./ CV;
    else
        % Alternative magnesium block (kept exactly; uses mldivide+warnings)
        mag_block = local_mag_block(x(i,:,1), exp(P.scale_NMDA));
        f(i,:,1) = ( GL*(VL - x(i,:,1)) + ...
                     x(i,:,2).*(VE - x(i,:,1)) + ...
                     x(i,:,3).*(VI - x(i,:,1)) + ...
                     x(i,:,5).*(VB - x(i,:,1)) + ...
                     x(i,:,6).*(VM - x(i,:,1)) + ...
                     x(i,:,7).*(VH - x(i,:,1)) + ...
                     x(i,:,4).*(VN - x(i,:,1)).*mag_block ) ./ CV;
    end
    

    % --- Conductance equations --------------------------------------------
    f(i,:,2) = (E'     - x(i,:,2)) .* (KE(i,:)');
    f(i,:,3) = (I'     - x(i,:,3)) .* (KI(i,:)');
    f(i,:,5) = (IB'    - x(i,:,5)) .* (KB(i,:)');
    f(i,:,4) = (ENMDA' - x(i,:,4)) .* (KN(i,:)');

    if IncludeMH
        f(i,:,6) = (Im' - x(i,:,6)) .* (KM(i,:) );
        f(i,:,7) = (Ih' - x(i,:,7)) .* (KH(i,:) );
    end

    % === Cross-talk across states (general operator) ======================
    % Flatten the *derivative* vector for this source
    f_i      = reshape(f(i,:,:), [np*nk, 1]);      % [56x1]
    df_xtalk = Xgain * (Kblock * f_i);             % cross-talk increment

    % Add cross-talk contribution back into the derivative
    f(i,:,:) = f(i,:,:) + reshape(df_xtalk, [1, np, nk]);
end

% --- Vectorise state derivatives ------------------------------------------
f = spm_vec(f);

% Pre-assign outputs for optional computation branches
[J,Q,D] = deal([]);

% --- Optional Jacobian -----------------------------------------------------
if (nargout < 2 || nargout == 50) && nargin < 5
    return
end

J = spm_cat(spm_diff(M.f, x, u, P, M, 1));

% --- Optional Delays -------------------------------------------------------
if nargout < 3 && nargin < 5
    return
end

% Fixed delay params (kept as in original logic)
D_   = [1 16];
d    = D_ .* full(exp(P.D(1:2))) / 1000;

% Same-source / same-population Kronecker scaffolds
Sp = kron(ones(nk,nk),kron( eye(np),eye(ns)));
Ss = kron(ones(nk,nk),kron(ones(np),eye(ns)));

% Cortico-thalamic / thalamo-cortical constants (ms -> s)
CT = 8;  % cortex->thalamus
TC = 3;  % thalamus->cortex

Tc              = zeros(np,np);
Tc([7 8],[1:6]) = CT * exp(P.CT); % L6->thal
Tc([1:6],[7 8]) = TC * exp(P.TC); % thal->cortex
Tc = Tc / 1000;
Tc = kron(ones(nk,nk),kron(Tc,ones(ns,ns)));

% Intra-population delays (ID), kept verbatim
ID = [2 1 1 1 1 2 1 2];
ID = ID .* exp(P.ID) / 1000;
ID = repmat(ID,[1 nk]);
ID = repmat(ID(:)',[np*nk,1]);
ID = kron(ID,ones(ns,ns));

% Historical alternatives kept off; final D as in original:
% D = d(1)*Ds + Tc + ID;  (overwritten in original)
D = Tc + ID;

% ---------------------- Local utilities (private) -------------------------
function mb = local_mag_block(vrow, scaleNMDA)
    % Preserves original mldivide + warning suppression semantics.
    % mb = 1 ./ (1 + 0.2*exp(-0.062*scaleNMDA * vrow));
    % Implemented via mldivide to match original exactly.
    wstate = warning('query','all'); %#ok<WNOFF>
    warning('off','all');
    denom  = (1 + 0.2*exp(-0.062*(scaleNMDA)*squeeze(vrow)'))';
    mb     = mldivide(denom, 1)';   % elementwise reciprocal via backslash
    [~,wid] = lastwarn;             %#ok<ASGLU>
    if ~isempty(wid), warning('off',wid); end
    warning(wstate); % restore
end

end
