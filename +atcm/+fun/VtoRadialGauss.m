function [G, b, GL, wid] = VtoRadialGauss(Q, w, model, ampwid)
    % Approximate conversion of vector Q to a symmetric radial feature matrix
    % with a smoother Gaussian kernel for each element in input Q.
    %
    % Usage: [G, b, GL, wid] = VtoRadialGauss(Q, w, model, ampwid)
    %
    % AS25
    
    if nargin < 4 || isempty(ampwid)
        ampwid = 1;
    end
    
    if nargin < 3 || isempty(model)
        model = 'Gauss';
    end
    
    if nargin < 2 || isempty(w)
        w = 4;
    end

    % Prerequisites
    nQ = length(Q);
    G = zeros(nQ, nQ);
    x = 1:nQ;
    Q = real(Q);

    % Scaled version of Q to estimate widths
    QQ = rescale(abs(Q), 1, (2.58) * 2);
    QQ = QQ .* sign(Q);
    wid = zeros(1, nQ);

    % Compute radial Gaussian kernel
    for i = 1:nQ
        % Central value and position
        v = Q(i);
        I = i;

        % Adjust width based on amplitude
        if ampwid && w == 4
            vv = QQ(i);
            w = max(real(vv), 1);    
            wid(i) = w;
        end

        % Pairwise radial Gaussian calculation
        for j = 1:nQ
            distance = abs(I - j); % Radial distance
            G(j, i) = v * exp(-distance^2 / (2 * w^2)); % Gaussian kernel
        end
    end

    % Make the matrix symmetric and normalize
    G = (G + G') / 2;
    G = G ./ norm(G, 'fro'); % Normalize by Frobenius norm

    % Rescale to min/max of input
    if min(Q) ~= max(Q)
        G = reshape(rescale(G(:), min(Q), max(Q)), size(G));
    end

    % Optional outputs
    if nargout == 2
        % Reduce if requested
        b = atcm.fun.lsqnonneg(G, diag(Q));
    else
        b = [];
    end

    if nargout >= 3
        % Graph Laplacian (optional output)
        Q = G;
        A = Q .* ~eye(length(Q));
        N = size(A, 1);
        GL = speye(N, N) + (A - spdiags(sum(A, 2), 0, N, N)) / 4;
    end
end
