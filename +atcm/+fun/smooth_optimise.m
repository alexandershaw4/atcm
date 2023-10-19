function y = smooth_optimise(y,yfix,tol)
% match the smoothness of vector y to that of vector yfix by optimising a
% smoothing/unsmoothing function that matches the norm of the jerkiness
% (third derivative). Def tol 1e-3.
%
%  dy = smooth_optimise(y,yfix,[tol]);
%
% AS2023

% compute jerk (~ approximate bumpyness)
g = @(x) norm(gradient(gradient(gradient(x))));

if nargin < 3 || isempty(tol); tol = 1e-3; end

iterate = true; n = 0;
while iterate
    if g(y) > g(yfix)
        % smooth it with a moving window
        y = atcm.fun.awinsmooth(y,4);
    elseif g(y) < g(yfix)
        % unsmooth it with its curvature
        y = y - 3*gradient(gradient(y));
    end

    % assess convergence
    if abs(g(yfix) - g(y)) < tol
        iterate = false;
    end

    % stop infinite runtime
    n = n + 1;
    if n == 1000
        iterate = false;
    end

end

end