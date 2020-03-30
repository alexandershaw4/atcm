function F = fe0(M,U,Y,PE,Gmax)
%
% log evidence
%--------------------------------------------------------------------------
% F   - [-ve] free energy F = log evidence = p(y|f,g,pE,pC) = p(y|m)
%
%__________________________________________________________________________

IS = M.IS;
y = Y.y;

% size of data (usually samples x channels)
%--------------------------------------------------------------------------
if iscell(y)
    ns = size(y{1},1);
else
    ns = size(y,1);
end
nr   = length(spm_vec(y))/ns;       % number of samples and responses
M.ns = ns;                          % store in M.ns for integrator

% initial states
%--------------------------------------------------------------------------
try
    M.x;
catch
    if ~isfield(M,'n'), M.n = 0;    end
    M.x = sparse(M.n,1);
end

% input
%--------------------------------------------------------------------------
try
    U;
catch
    U = [];
end

% initial parameters
%--------------------------------------------------------------------------
M.P = M.pE;

% time-step
%--------------------------------------------------------------------------
try
    Y.dt;
catch
    Y.dt = 1;
end

% precision components Q
%--------------------------------------------------------------------------
try
    Q = Y.Q;
    if isnumeric(Q), Q = {Q}; end
catch
    Q = spm_Ce(ns*ones(1,nr));
end
nh    = length(Q);                  % number of precision components
nt    = length(Q{1});               % number of time bins
nq    = nr*ns/nt;                   % for compact Kronecker form of M-step
h     = zeros(nh,1);                % initialise hyperparameters


% confounds (if specified)
%--------------------------------------------------------------------------
try
    nb   = size(Y.X0,1);            % number of bins
    nx   = nr*ns/nb;                % number of blocks
    dfdu = kron(speye(nx,nx),Y.X0);
catch
    dfdu = sparse(ns*nr,0);
end

% hyperpriors - expectation
%--------------------------------------------------------------------------
try
    hE = M.hE;
    if length(hE) ~= nh
        hE = hE*sparse(nh,1);
    end
catch
    hE = sparse(nh,1);
end

% prior moments
%--------------------------------------------------------------------------
pE    = M.pE;
pC    = M.pC;
nu    = size(dfdu,2);                 % number of parameters (confounds)
np    = size(pC,2);                   % number of parameters (effective)

% [Alex]: unpack potentially structural (co)variance
%---------------------------------------------------
if isstruct(pC)
    pC = diag(spm_vec(pC));
end

% second-order moments (in reduced space)
%--------------------------------------------------------------------------
ipC   = spm_inv(pC);

% initialize conditional density
%--------------------------------------------------------------------------
Eu    = spm_pinv(dfdu)*spm_vec(y);
Ep    = pE;

% precision and conditional covariance
%------------------------------------------------------------------
iS    = sparse(0);
for i = 1:nh
    iS = iS + Q{i}*(exp(-16) + exp(hE(i)));
end
S     = spm_inv(iS);
iS    = kron(speye(nq),iS);
qE    = spm_vec(pE);
pE    = spm_vec(pE);
y     = spm_vec(y);
np    = length(qE);

% Sampling
%==========================================================================
if nargin < 5 || isempty(Gmax)
    Gmax  = -Inf;
end

% new param set prediction
%------------------------------------------------------------------
R = spm_vec(feval(IS,spm_unvec(spm_vec(PE),M.pE),M,U));
P = spm_vec(PE);

% prediction error
%------------------------------------------------------------------
ey     = R - y;  % spectral (objective) error
ep     = P - pE; % parameter error
ep(isinf(ep)) = 0;

% Gibb's energy
%------------------------------------------------------------------
qh     = real(ey')*iS*real(ey) + imag(ey)'*iS*imag(ey);
G(i,1) = - ns*log(qh)/2 - ep'*ipC*ep/2;

% conditional mode
%----------------------------------------------------------------------
[maxG,j] = max(G);
if maxG  > Gmax
    qE   = P(:,j);
    f    = R(:,j);
    Gmax = maxG;
end
    
% conditional dispersion
%----------------------------------------------------------------------
q     = exp((G - maxG));
q     = q/sum(q);
for i = 1:np
    for j = 1:np
        qC(i,j) = ((P(i,:) - qE(i)).*(P(j,:) - qE(j)))*q;
    end
end


% objective function:
%======================================================================
F     = Gmax + spm_logdet(ipC*qC)/2;
F     = Gmax;


    
