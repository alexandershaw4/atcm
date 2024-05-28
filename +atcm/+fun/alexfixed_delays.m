function x0 = alexfixed_delays(P,M,tol,a,input,dt)
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

if nargin < 6 || dt == 0
    dt = 1;
end

if nargin < 5 || isempty(input)
    input = 0;
end

f  = @(x,varargin) M.f(x,input,P,M);
x0 = M.x(:);

if nargin < 3 || isempty(tol)
    tol = 1e-6;%10;
end

if nargin < 4 || isempty(a)
    a = 2;
end

% fetch delay operator
[dx,A,D] = f(x0,[]);
D        = real(D);

for i = 1:3e3

    [dx,A] = f(x0,[]);

    % delay model
    
    % work in seconds;
    dx = dx / dt;

    % compute a linear update operator for this step
    b    = pinv(full(A)'.*x0).*dx;
    Q    = (A.*b); % dxdt = Q*x;

    % add delays and recompute step @ dt
    Q    = Q + (D.*~~real(A));;
    dx   = (dt*Q*x0);

    x1   = x0 - (exp(-a) * pinv(full(A)) * dx);
    e(i) = norm(f(x1,[]) );
    
    norm(f(x1,[]) )

    if norm(f(x1,[]) ) < tol
        %fprintf('Converged at iteration %d\n', i);
        break;
    end

    if i > 1 && e(i) > e(i-1)
        a = a * 2;
    end

    x0 = x1;

end

