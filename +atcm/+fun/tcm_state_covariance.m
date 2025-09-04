function [Sigma_x, J, info] = tcm_state_covariance(tcm_fun, xstar, ustar, P, M, Q)
% Compute steady-state state covariance from the local Jacobian of a TCM.
%
% Inputs
%   tcm_fun : handle to your TCM state-eq function, e.g. @tcm2025a
%             must return [f,J,D] with J = df/dx evaluated at (x,u,P,M)
%   xstar   : fixed point (or operating point) state vector (ns x 1)
%   ustar   : input vector at the operating point (nu x 1), can be zeros(ns,1)
%   P, M    : parameter and model structs for your TCM
%   Q       : process noise covariance (ns x ns). If scalar, uses Q*I.
%
% Outputs
%   Sigma_x : steady-state state covariance solving  J*Σ + Σ*J' + Q = 0
%   J       : Jacobian at (xstar, ustar)
%   info    : struct with diagnostics (eigenvalues, stable flag, msg)
%
% Notes
% - Requires a stable linearisation (Re(lambda(J)) < 0) for Lyapunov solution.
% - For delay systems, this treats J as the no-delay linearisation.
%   (Use your Laplace-domain transfer if you want delay-aware spectra.)
%
% Example:
%   [f0,J] = tcm2025a(xstar, ustar, P, M);
%   Q = 1e-4 * eye(numel(xstar));
%   [Sigma_x, J] = tcm_state_covariance(@tcm2025a, xstar, ustar, P, M, Q);

    if isscalar(Q), Q = Q * eye(numel(xstar)); end

    % Evaluate model and Jacobian at operating point
    [~, J] = tcm_fun(xstar, ustar, P, M);
    J = full(J);

    % Diagnostics
    lam = eig(J);
    stable = all(real(lam) < 0);
    info = struct('eigJ', lam, 'stable', stable, ...
                  'msg', ternary(stable, 'OK: stable linearisation', ...
                                         'WARNING: unstable linearisation'));

    % Solve continuous-time Lyapunov: J*Σ + Σ*J' + Q = 0
    % MATLAB lyap solves A*X + X*A' + Q = 0
    Sigma_x = lyap(J, Q);

    labels = {'V' 'AMPA' 'GABAA' 'NMDA' 'GABAB' 'M' 'H'};
    groups = {1:8, 9:16, 17:24, 25:32, 33:40, 41:48, 49:56};

    plot_state_covariance(Sigma_x, groups, labels)

end

function y = ternary(cond, a, b)
    if cond, y = a; else, y = b; end
end

function plot_state_covariance(Sigma, groups, labels)
% Plot covariance matrix with group boxes.
%
% Inputs
%   Sigma  : covariance matrix (n x n)
%   groups : cell array, each element is vector of indices for that group
%   labels : optional cell array of group labels (for axes)

    figure;
    imagesc(Sigma);
    axis square; colormap(parula); colorbar;
    hold on;

    n = size(Sigma,1);

    % draw boxes for each group
    for g = 1:numel(groups)
        idx = groups{g};
        i1 = min(idx); 
        i2 = max(idx);

        % draw rectangle: [x,y,w,h]
        rectangle('Position',[i1-0.5, i1-0.5, i2-i1+1, i2-i1+1], ...
                  'EdgeColor','k','LineWidth',1.5);
    end

    % optional: add group separator lines
    for g = 1:numel(groups)
        i2 = max(groups{g});
        xline(i2+0.5,'k-','LineWidth',1);
        yline(i2+0.5,'k-','LineWidth',1);
    end

    % optional: label ticks
    if nargin > 2
        centers = cellfun(@(ix) mean([min(ix) max(ix)]), groups);
        set(gca,'XTick',centers,'YTick',centers,'XTickLabel',labels,'YTickLabel',labels);
        xtickangle(45);
    end

    title('State covariance with group boxes');
end
