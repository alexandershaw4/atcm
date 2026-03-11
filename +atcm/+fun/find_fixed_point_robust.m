function x_fp = find_fixed_point_robust(P, M, tol, max_iter, verbose, input)
% Robust fixed-point solver for tc_hilge2 / DCM RHS using
% damped Gauss-Newton + backtracking line-search on ||f||.
% Falls back to a quasi-Newton root-finder if needed.
%
% Inputs
%   P, M     : model structs
%   x0       : initial state (struct like M.x OR vector)
%   tol      : tolerance on ||f||
%   max_iter : max iterations (default 200)
%   verbose  : 0/1
%   input    : scalar or vector input (default 0)
%
% Output
%   x_fp     : fixed point (same shape as M.x if possible)

if nargin < 4 || isempty(tol),      tol = 1e-6; end
if nargin < 5 || isempty(max_iter), max_iter = 200; end
if nargin < 6 || isempty(verbose),  verbose = true; end
if nargin < 7 || isempty(input),    input = 0; end

x0 = M.x(:);

% --- vectorisation helpers (SPM-style if available)
use_spm = exist('spm_vec','file') == 2 && isfield(M,'x');

if use_spm
    vec   = @(xs) spm_vec(xs);
    unvec = @(xv) spm_unvec(xv, M.x);
else
    % fallback: assume x0 is already a vector
    vec   = @(xs) xs(:);
    unvec = @(xv) xv(:);
end

% --- evaluate f and J at a vector state
    function [fvec, J] = fun_and_jac(xv)
        xs = unvec(xv);

        % Try to call a tc_hilge2-style function if M.f is not provided
        if isfield(M,'f') && ~isempty(M.f)
            % Expect [f, J] = M.f(x,u,P,M) (SPM/DCM convention)
            [f_out, J_out] = M.f(xs, input, P, M);
        else
            % If you want this specifically for tc_hilge2, replace this block:
            error('M.f not found; please provide M.f that returns [f,J].');
        end

        fvec = vec(f_out);
        J    = full(J_out); % ensure numeric full matrix
    end

% --- init
x = vec(x0);
[f, J] = fun_and_jac(x);
nf = norm(f);

% variable scaling: scale_i = 1/(|x_i|+1), clipped
scale = 1 ./ (abs(x) + 1);
scale = max(min(scale, 1e3), 1e-3);

if verbose
    fprintf('[fp] init ||f|| = %.3e\n', nf);
end

% --- main loop: damped Gauss-Newton on ||f||
for it = 0:(max_iter-1)

    if ~isfinite(nf)
        if verbose, fprintf('[fp] non-finite ||f||; bailing to fallback.\n'); end
        break;
    end

    if nf < tol
        if verbose
            fprintf('[fp] converged at iter %d, ||f||=%.3e\n', it, nf);
        end
        x_fp = unvec(x);
        return;
    end

    % Column scaling: Js = J ./ scale'  (scale columns)
    Js = J ./ (scale(:)');   % size: [nEq x nVar]
    fs = f;

    % Ridge similar to python:
    lam = 1e-8 + 1e-6 * (nf / (norm(x) + 1e-12));

    % Solve (Js'Js + lam I) dxs = -Js' f
    A = (Js.' * Js) + lam * speye(size(Js,2));
    b = -(Js.' * fs);

    % Use backslash (sparse-safe); convert to full if needed
    dxs = A \ b;

    % Unscale
    dx = dxs ./ scale(:);

    % Backtracking line search on ||f||
    t = 1.0;
    nf0 = nf;
    accept = false;

    for ls = 1:25
        x_new = x + t * dx;
        [f_new, J_new] = fun_and_jac(x_new);
        nf_new = norm(f_new);

        % Armijo-ish decrease on norm (same as python logic)
        if isfinite(nf_new) && (nf_new <= nf0 * (1 - 1e-4 * t) || nf_new < nf0)
            x = x_new;
            f = f_new;
            J = J_new;
            nf = nf_new;
            accept = true;
            break;
        end

        t = t * 0.5;
    end

    if verbose
        if accept, msg = 'ok'; else, msg = 'fail'; end
        fprintf('[fp] it %03d ||f||=%.3e -> %.3e  step=%.2e  %s\n', it, nf0, nf, t, msg);
    end

    if ~accept
        % Newton couldn't find a decreasing step → fallback
        break;
    end

    % small-step termination
    if norm(t * dx) < 1e-10 * (norm(x) + 1.0)
        if verbose, fprintf('[fp] step too small; bailing to fallback.\n'); end
        break;
    end
end

% --- fallback: quasi-Newton root solve on f(x)=0 (no need for J)
% We minimise 0.5||f||^2 with fminunc if Optimization TB exists; otherwise use fminsearch.

fun_only = @(xv) vec(fun_and_jac(xv)); % WARNING: this returns [f; J]? so don't use
% define a clean fun-only:
fun_only = @(xv) local_fun_only(xv);

    function fvec = local_fun_only(xv)
        [fvec, ~] = fun_and_jac(xv);
    end

phi = @(xv) 0.5 * (local_fun_only(xv).' * local_fun_only(xv));

if verbose
    fprintf('[fp] entering fallback...\n');
end

x_try = x;

if license('test','Optimization_Toolbox') && exist('fminunc','file') == 2
    opts = optimoptions('fminunc', 'Display','off', 'MaxIterations',200, ...
        'OptimalityTolerance', tol, 'StepTolerance', 1e-12);
    try
        x_try = fminunc(phi, x_try, opts);
    catch
        % ignore
    end
else
    opts = optimset('Display','off', 'MaxIter',200, 'TolX',1e-12, 'TolFun',tol);
    try
        x_try = fminsearch(phi, x_try, opts);
    catch
        % ignore
    end
end

res_norm = norm(local_fun_only(x_try));

if verbose
    fprintf('[fp] final residual ||f|| = %.3e\n', res_norm);
end

x_fp = unvec(x_try);

% safety checks (same spirit as python)
bad = (~isfinite(res_norm)) || (res_norm > max(1e-2, 10*tol)) || ...
      any(~isfinite(vec(x_fp))) || (norm(vec(x_fp)) > 1e8);

if bad
    if verbose, fprintf('[fp] rejecting fixed point; returning x0.\n'); end
    x_fp = x0;
end

end
