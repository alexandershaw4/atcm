function y = aINTw_ica(P,M,U)

pE  = P;                          % real world parameter vector
iW  = M.iW;                       % mapping from subspace 
C   = M.C;                        % component vector

x1 = (pE*iW) .* spm_vec(M.p0E)';

% x0  = (diag(pE) * hpE);           % indexing
% x0  = sum(x0);
% x0(x0==0) = 1;

%x1  = x0.*spm_vec(M.p0E)';   % new (full) parameter values
PP  = spm_unvec(x1,M.p0E);        % proper structure

M.pC = ( ~~spm_vec(PP) +0.125 );
M.pC = spm_unvec(M.pC,PP);

IS = spm_funcheck(M.f0);       % Integrator
y  = IS(PP,M,U);                  % Prediction

