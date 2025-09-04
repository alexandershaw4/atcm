function [Syy, Sxx] = tcm_psd_from_jacobian(J, Q, C, R, freqs)
% Compute theoretical PSD from local linearisation (ignoring delays).
% freqs in Hz; returns Sxx(:,:,k) and Syy(:,:,k) for each frequency.
    if isscalar(Q), Q = Q * eye(size(J,1)); end
    if isscalar(R), R = R * eye(size(C,1)); end

    ns = size(J,1); ny = size(C,1);
    K  = numel(freqs);
    Sxx = zeros(ns, ns, K);
    Syy = zeros(ny, ny, K);

    for k = 1:K
        w  = 2*pi*freqs(k);
        H  = ((1i*w)*eye(ns) - J) \ eye(ns);  % (iωI - J)^{-1}
        Sx = H * Q * (H');                    % H Q H^H
        Sxx(:,:,k) = real(Sx);                % PSD is Hermitian; keep real part
        Syy(:,:,k) = real(C * Sx * C' + R);
    end
end
