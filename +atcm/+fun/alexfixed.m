function x0 = alexfixed(P,M,tol,a,input,numit)
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

if nargin < 5 || isempty(numit)
    numit = 3e3;    
end

if nargin < 5 || isempty(input)
    input = 0;
end

%if dt > 0
 %   f  = @(x,varargin) M.f(x,input,P,M,dt);
%else
    f  = @(x,varargin) M.f(x,input,P,M);
%end

x0 = M.x(:);

if nargin < 3 || isempty(tol)
    tol = 1e-6;%10;
end

if nargin < 4 || isempty(a)
    a = 2;
end

% fetch delay operator
%[dx,A,D] = f(x0,[]);
%D        = inv(eye(length(D)) - D);

%rf = 0.2; 

for i = 1:numit;%3e10

    [dx,A] = f(x0,[]);

    %A = D*A;
    
    x1 = x0 - (exp(-a) * pinv(full(A)) * dx);

    % x1 = x0 - dt * pinv(full(A)) * dx;

    %x1 = ( (1 - rf) * x0 ) - ( rf * exp(-a) * pinv(full(A)) * dx );

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

