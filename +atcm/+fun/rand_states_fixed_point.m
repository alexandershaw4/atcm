function [xfp,dx] = rand_states_fixed_point(P,M)
% Robust Fixed Point
% iteratively find a fixed point in the model under randomised initial
% conditions
% [xfp,dx] = atcm.fun.rand_states_fixed_point(P,M)
% AS22

x = M.x;
v = abs(x/64);
v(v==0) = 1/64;
k = 1e-12;

if any(M.u) && length(M.u) > 1
    %M.u = min(M.u);
    M.u = min(M.u(M.u~=0));
end

for i = 1:24
    
    % Laplacian sampling
    r   = ( x + (v).*rand(size(x)) );

    M.x = r;
        
    % Solve
    x   = solvefixed(P,M,k);
    
    % Store
    dx(:,i) = x(:);
    
end

% Average all fps
dx(isnan(dx))=0;
xfp = spm_unvec(mean((dx),2),M.x);

end


function [x] = solvefixed(P,M,k)

% solve for fixed point
%------------------------------------------------------------------
ns    = size(P.A{1},1);     % number of sources (endogenous inputs)
a     = 2;                  % regulariser
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
    if dnx < k; break, end 
    
end
x = M.x;
end