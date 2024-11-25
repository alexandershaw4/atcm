function [X,F] = agd(f,x,n,a)
% super simple gradient descent of a multivariate function using Newton
% steps, where the step size is controlled by the curvature (Newton's
% method).
%
% [X,F] = agd(fun,x0,Niterations)
%
% AS

e = f(x);

if nargin < 4 || isempty(a);
    a = 1;
end

doplot = 1;
alle   = e;

fprintf('Initialise\n');
fprintf('It: %d | f = %d\n',0,e);

% iterate
for i = 1:n
    
    % gradients
    g = dfdx(f,x);

    % check whether f'(x) == f(x)
    if norm(g) < 1e-6
        fprintf('finished\n');
        X = x; F = e;
        return;
    end

    % if n == 1, optimise alpha
    %if i == 1
    %    fa = @(a) f(x + (a)*-g(:));
    %    ga = dfdx(fa,a);
    %    a  = e./ga;
    %end

    % if i == 1
    %     a = 1/8;
    % end

    H = dfdxdx(f,x);

    %H = g(:)*g(:)';

    g = g./norm(g);
    H = H./norm(H);
    H = makeposdef(H);

    % step and evaluate
    x  = x - a*(H\g');
    
    de = f(x);

    if de > e
        a = a ./ 2;
    end

    e = de;

    fprintf('It: %d | f = %d\n',i,e);
    alle = [alle e];

    % stop at maximum iterations
    if i == n
        X = x;
        F = e;
        fprintf('Did not finish\n');
    end

    if doplot
        plot(1:i+1,alle,'*',1:i+1,alle);drawnow;
    end

end


end

function g = dfdxdx(f,x)
% simple 1-step finite difference routine for compute partial gradients

e0 = f(x);
k  = exp(-8);

for i  = 1:length(x)
    for j = 1:length(x)
        if i < j
            dx    = x;
            dx(i) = dx(i) + k;
            dx(j) = dx(j) + k;
        
            g(i,j)  = (f(dx) - e0) ./ (2*k);
            g(j,i)  = g(i,j);
        end
    end
end

end

function g = dfdx(f,x)
% simple 1-step finite difference routine for compute partial gradients

e0 = f(x);
k  = exp(-8);

for i  = 1:length(x)
    dx    = x;
    dx(i) = dx(i) + k;
    g(i)  = (f(dx) - e0) ./ k;
end

end