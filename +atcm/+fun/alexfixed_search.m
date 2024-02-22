function x0 = alexfixed_search(P,M,tol,a,input)
% Use fminsearch to find a fixed point of a dynamical system
% (in this case, a DCM) with fixed step size (regulariser);
%
% x0 = atcm.fun.alexfixed_search(P,M)
%
% P = (DCM) parameter structure
% M = (DCM) model structure with;
% - M.f = handle to rhs (model function)
% - M.x = initial states such that the dynamical model is specified:
%    [dx,dfdx] = M.f(M.x,0,P,M);
%
% see alexfixed for the Newton Raphson algorithm version
%
% AS

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


[x0,FVAL]= fminsearch(@(x) fun(f,x),[-2; x0]);

end

function e = fun(f,x0)

step = x0(1);
x0 = x0(2:end);

[dx,A] = f(x0,[]);

x1 = x0 - exp(step) * pinv(full(A)) * dx;

e = sum( (x1 - x0).^2 );

end


% for i = 1:3e10
% 
%     [dx,A] = f(x0,[]);
% 
% 
%     x1 = x0 - exp(-a) * pinv(full(A)) * dx;
% 
%     %x1 = ( (1 - rf) * x0 ) - ( rf * exp(-a) * pinv(full(A)) * dx );
% 
%     e(i) = norm(x1 - x0);
% 
%     norm(x1 - x0)
% 
%     if norm(x1 - x0) < tol
%         %fprintf('Converged at iteration %d\n', i);
%         break;
%     end
% 
%     if i > 1 && e(i) > e(i-1)
%         a = a * 2;
%     end
% 
%     x0 = x1;
% 
% end

