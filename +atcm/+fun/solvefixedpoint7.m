function [x,f] = solvefixedpoint7(P,M)

% inital states
%------------------------------------------------------------------
ns   = size(P.A{1},1);                           % number of sources
np   = 6;                                        % number of populations

% create (initialise voltage at -50mV)
%--------------------------------------------------------------------------
x        = zeros(ns,np,5) ;
x(:,:,1) = -70;
        
M.g   = {};
M.x   = x;
M.pE  = P;
M.n   = length(spm_vec(x));
M.u   = sparse(ns,1);

% solve for steady state
%--------------------------------------------------------------------------
warning off
x     = solvefixed(P,M);
f     = M.f;
warning on
end

function x = solvefixed(P,M)

% solve for fixed point
%------------------------------------------------------------------
ns    = size(P.A{1},1);     % number of sources (endogenous inputs)
a     = 2;                  % regulariser
M.u   = sparse(ns,3);
dnx   = 0;
for i = 1:128
    
    % solve under locally linear assumptions
    %--------------------------------------------------------------
    [f,dfdx] = feval(M.f,M.x,M.u,P,M);
    dx       = - dfdx\f;
    
    % regularise
    %--------------------------------------------------------------
    ndx   = norm(dx,Inf);
    if ndx < dnx
        a = a/2;
    end
    dnx    = ndx;
    
    % update and convergence
    %--------------------------------------------------------------
    M.x    = spm_unvec(spm_vec(M.x) + exp(-a)*dx,M.x);
    if dnx < 1e-12, break, end
    
end

x = M.x;
end