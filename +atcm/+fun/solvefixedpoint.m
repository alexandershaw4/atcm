function [x,i] = solvefixedpoint(P,M,input)


% inital states
%------------------------------------------------------------------
ns   = size(P.A{1},1);                           % number of sources

if isfield(M,'x') && ~iscell(M.x) && ~isempty(M.x)
    np = size(M.x,2);
else
    np   = size(P.H,1);                              % number of populations
end

%np = 7;

% find fp in presence of an input (DC), or not
if nargin < 3
    M.u = sparse(ns,1);%3
else
    M.u = input(1);
end
    
% create (initialise voltage at -50mV)
%--------------------------------------------------------------------------
if isfield(M,'x') && ~iscell(M.x)
    nk = size(M.x,3);
else
    nk = 7;
end

%if isfield(M,'x') && ~iscell(M.x)
    % Point estimator / neural mass models
    mfm = 0;
    x        = zeros(ns,np,nk) ;
    x(:,:,1) = -70;
%elseif isfield(M,'x') && iscell(M.x)
%    % Mean field models
%    mfm = 1;
%    x{1}        = zeros(ns,np,nk) ;
%    x{1}(:,:,1) = -70;
%    x{2}        = ones(nk,nk,ns,np);
%end
    
M.g   = {};
M.x   = x;
M.pE  = P;
M.n   = length(spm_vec(x));

% solve for steady state
%--------------------------------------------------------------------------
warning off
[x,i] = solvefixed(P,M,mfm);
f     = M.f;
warning on

end

function [x,i] = solvefixed(P,M,mfm)

% solve for fixed point
%------------------------------------------------------------------
ns    = size(P.A{1},1);     % number of sources (endogenous inputs)
a     = 2;                  % regulariser
dnx   = 0;
for i = 1:128
    
    % solve under locally linear assumptions
    %--------------------------------------------------------------
    [f,dfdx] = feval(M.f,M.x,M.u,P,M);
    
%     if mfm
%         mf = spm_unvec(f,M.x);
%         f  = spm_vec(mf{1});
%     end
    
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
    if dnx < 1e-6, break, end % 12
    
end
x = M.x;
end