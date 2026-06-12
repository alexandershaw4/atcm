function [f, J, D, aux] = tc_ngnmm_fx(x, u, P, M)
% tc_ngnmm_fx
% 8-population next-generation neural mass model in DCM-compatible format.
%
% Accepted input state formats:
%   x = Ns x Npop x 4
%   x = spm_vec(M.x)-ordered vector, which is reshaped using M.x
%
% State dimension 3:
%   1 = R : population firing rate
%   2 = V : mean voltage
%   3 = g : lumped incoming conductance
%   4 = h : conductance derivative
%
% Outputs:
%   f   = vectorised dx/dt, using x(:)/spm_vec-compatible ordering
%   J   = d vec(f) / d vec(x), analytic state Jacobian
%   D   = state-by-state delay matrix in seconds
%   aux = optional diagnostics

if nargin < 2 || isempty(u)
    u = 0;
end

% ---------------------------------------------------------------------
% Reshape vector input back to DCM.M.x shape if needed
% ---------------------------------------------------------------------
if isvector(x)
    if isfield(M,'x') && numel(x) == numel(M.x)
        x = reshape(x, size(M.x));
    else
        error('tc_ngnmm_fx:BadState', ...
            'Vector x has %d elements but cannot be reshaped using M.x.', numel(x));
    end
end

sx = size(x);
if numel(sx) < 3 || sx(3) ~= 4
    error('tc_ngnmm_fx expects x to be Ns x Npop x 4, or a vector matching M.x.');
end

Ns = sx(1);
Np = sx(2);
Nx = sx(3);
nvec = numel(x);

if isfield(M,'Npop') && M.Npop ~= Np
    error('tc_ngnmm_fx:MismatchedNpop', ...
        'M.Npop=%d but size(x,2)=%d.', M.Npop, Np);
end

% ---------------------------------------------------------------------
% Indices into vectorised state
% ---------------------------------------------------------------------
idx = @(s,p,k) sub2ind(sx,s,p,k);

% ---------------------------------------------------------------------
% Robust parameter defaults
% ---------------------------------------------------------------------
if ~isfield(P,'logDelta'), P.logDelta = log(0.5) * ones(1,Np); end
if ~isfield(P,'eta'),      P.eta      = zeros(1,Np);          end
if ~isfield(P,'logAlpha'), P.logAlpha = log(100) * ones(1,Np);end

if ~isfield(P,'logKappa')
    if isfield(M,'Kmask')
        tmp = log(1e-6 * ones(Np,Np));
        tmp(M.Kmask == 1) = log(0.2);
        P.logKappa = tmp;
    else
        P.logKappa = log(0.2 * ones(Np,Np));
    end
end

if ~isfield(P,'Erev')
    if isfield(M,'Erev')
        P.Erev = M.Erev;
    else
        P.Erev = 60 * ones(Np,Np);
        if isfield(M,'popType')
            P.Erev(:, M.popType < 0) = -90;
        elseif Np >= 3
            P.Erev(:, [3 5 7]) = -90;
        end
    end
end

if ~isfield(P,'inputGain'), P.inputGain = zeros(1,Np); end

if isfield(P,'logQIFRate')
    qifRate = exp(P.logQIFRate);
else
    qifRate = 1;
end

Delta = exp(P.logDelta(:))';
eta   = P.eta(:)';
alpha = exp(P.logAlpha(:))';
Bin   = P.inputGain(:)';

if numel(Delta) ~= Np, Delta = Delta(1) * ones(1,Np); end
if numel(eta)   ~= Np, eta   = eta(1)   * ones(1,Np); end
if numel(alpha) ~= Np, alpha = alpha(1) * ones(1,Np); end
if numel(Bin)   ~= Np, Bin   = zeros(1,Np);           end

Kappa = exp(P.logKappa);
if isvector(Kappa)
    if numel(Kappa) == Np*Np
        Kappa = reshape(Kappa(:), Np, Np);
    else
        Kappa = Kappa(1) * ones(Np,Np);
    end
end

if isfield(M,'Kmask')
    Kappa = Kappa .* M.Kmask;
end

Erev = P.Erev;
if isvector(Erev)
    if numel(Erev) == Np*Np
        Erev = reshape(Erev(:), Np, Np);
    else
        Erev = Erev(1) * ones(Np,Np);
    end
end

% Effective reversal per target population, weighted by current gains.
Eeff = zeros(1,Np);
for p = 1:Np
    denom = sum(Kappa(p,:)) + eps;
    Eeff(p) = sum(Kappa(p,:) .* Erev(p,:)) ./ denom;
end

% ---------------------------------------------------------------------
% Unpack states
% ---------------------------------------------------------------------
R = x(:,:,1);
V = x(:,:,2);
g = x(:,:,3);
h = x(:,:,4);

Rraw = R;
R = max(R, 1e-8);

% External input
Iext = zeros(Ns,Np);
if ~isempty(u)
    try
        uu = full(u(1));
    catch
        uu = 0;
    end
    Iext = repmat(Bin, Ns, 1) .* uu;
end

% ---------------------------------------------------------------------
% Intrinsic target x source conductance drive within each source/region
% ---------------------------------------------------------------------
drive_g = zeros(Ns,Np);
for s = 1:Ns
    % Kappa rows = target populations, columns = source populations.
    drive_g(s,:) = (Kappa * R(s,:)')';
end

% ---------------------------------------------------------------------
% NGNMM dynamics
% ---------------------------------------------------------------------
dR = zeros(Ns,Np);
dV = zeros(Ns,Np);
dg = h;
dh = zeros(Ns,Np);

for s = 1:Ns
    for p = 1:Np
        dR(s,p) = qifRate * ( ...
            Delta(p)/pi ...
            + 2 * R(s,p) * V(s,p) ...
            - g(s,p) * R(s,p) );

        dV(s,p) = qifRate * ( ...
            V(s,p)^2 ...
            + eta(p) ...
            + Iext(s,p) ...
            - (pi * R(s,p))^2 ...
            + g(s,p) * Eeff(p) ...
            - g(s,p) * V(s,p) );

        dh(s,p) = alpha(p)^2 * (drive_g(s,p) - g(s,p)) ...
            - 2 * alpha(p) * h(s,p);
    end
end

f3 = zeros(sx);
f3(:,:,1) = dR;
f3(:,:,2) = dV;
f3(:,:,3) = dg;
f3(:,:,4) = dh;

% Return vectorised derivative for compatibility with transfer functions.
f = f3(:);

% ---------------------------------------------------------------------
% Analytic Jacobian J = d vec(f) / d vec(x)
% ---------------------------------------------------------------------
J = zeros(nvec,nvec);

for s = 1:Ns
    for p = 1:Np
        iR = idx(s,p,1);
        iV = idx(s,p,2);
        ig = idx(s,p,3);
        ih = idx(s,p,4);

        % dR = Delta/pi + 2RV - gR
        J(iR,iR) = 2 * V(s,p) - g(s,p);
        J(iR,iV) = 2 * R(s,p);
        J(iR,ig) = -R(s,p);

        % dV = V^2 + eta + I - (pi R)^2 + gE - gV
        J(iV,iR) = -2 * pi^2 * R(s,p);
        J(iV,iV) =  2 * V(s,p) - g(s,p);
        J(iV,ig) =  Eeff(p) - V(s,p);

        % dg = h
        J(ig,ih) = 1;

        % dh = alpha^2 * (sum_q Kappa(p,q) R_q - g) - 2 alpha h
        for q = 1:Np
            iRq = idx(s,q,1);
            J(ih,iRq) = J(ih,iRq) + alpha(p)^2 * Kappa(p,q);
        end
        J(ih,ig) = -alpha(p)^2;
        J(ih,ih) = -2 * alpha(p);
    end
end

% If R was clipped, zero derivatives wrt clipped firing-rate states.
clipped = Rraw <= 1e-8;
if any(clipped(:))
    [ss, pp] = find(clipped);
    for n = 1:numel(ss)
        J(:, idx(ss(n),pp(n),1)) = 0;
    end
end

% Scale QIF rows of the Jacobian to match dR/dV scaling
for s = 1:Ns
    for p = 1:Np
        iR = idx(s,p,1);
        iV = idx(s,p,2);

        J(iR,:) = qifRate * J(iR,:);
        J(iV,:) = qifRate * J(iV,:);
    end
end

% ---------------------------------------------------------------------
% State-by-state delay matrix D
% ---------------------------------------------------------------------
D = zeros(nvec,nvec);

if isfield(P,'logDelay')
    baseDelay = exp(P.logDelay);
elseif isfield(P,'d')
    baseDelay = exp(P.d);
elseif isfield(M,'delay')
    baseDelay = M.delay;
else
    baseDelay = 0;
end
if numel(baseDelay) > 1, baseDelay = baseDelay(1); end

% Delays are placed on R_pre -> h_post dependencies. This is the
% delayed path by which pre-synaptic firing drives post-synaptic alpha
% conductance dynamics.
for s = 1:Ns
    for p = 1:Np
        ih = idx(s,p,4);
        for q = 1:Np
            if Kappa(p,q) ~= 0
                iRq = idx(s,q,1);
                D(ih,iRq) = baseDelay;
            end
        end
    end
end

% ---------------------------------------------------------------------
% Optional diagnostics
% ---------------------------------------------------------------------
if nargout > 3
    W = pi .* R + 1i .* V;
    Z = (1 - conj(W)) ./ (1 + conj(W));

    aux = struct();
    aux.f3 = f3;
    aux.R = R;
    aux.V = V;
    aux.g = g;
    aux.h = h;
    aux.W = W;
    aux.Z = Z;
    aux.sync = abs(Z);
    aux.phase = angle(Z);
    aux.Kappa = Kappa;
    aux.Erev = Erev;
    aux.Eeff = Eeff;
    aux.drive_g = drive_g;
    aux.D = D;
end
end
