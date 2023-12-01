function x0 = alexfixed(P,M)
% Use the Newton-Raphson method to find a fixed point of a dynamical system
% (in this case, a DCM) with adaptive step size (regulariser) that
% increases when iterations don't improve
%
% x0 = alexfixed(P,M)
%
% P = (DCM) parameter structure
% M = (DCM) model structure with;
% - M.f = handle to rhs (model function)
% - M.x = initial states such that the dynamical model is specified:
%    [dx,dfdx] = M.f(M.x,0,P,M);
%
% AS


f  = @(x,varargin) M.f(x,0,P,M);
x0 = M.x(:);
a  = 2;

tol = 1e-6;

for i = 1:1e10

    [dx,A] = f(x0,[]);
    
    x1 = x0 - exp(-a) * pinv(full(A)) * dx;

    e(i) = norm(x1 - x0);

    %norm(x1 - x0)

    if norm(x1 - x0) < tol
        fprintf('Converged at iteration %d\n', i);
        break;
    end

    if i > 1 && e(i) > e(i-1)
        a = a * 2;
    end

    x0 = x1;

end

