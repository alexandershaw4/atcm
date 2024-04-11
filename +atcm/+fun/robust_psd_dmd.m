function Ppf = robust_psd_dmd(x,dt,w,N)
% compute the PSD of the signal x using DMD of the Hilbert time-frequency
% representation; returns PSd at freqs in w with roughly N-components.
%
% Pf = atcm.fun.robust_psd_dmd(signal,dt,w,N)
%
% AS2024

[~,HTF,TFD] = atcm.fun.tfdecomp(x,dt,w);

X=TFD(:,1:end-1);
Y=TFD(:,2:end);

[U, S, V]=svd(X,'econ');

U=U(:,1:N);
V=V(:,1:N);
S=S(1:N,1:N);

A_tilde=U'*Y*V/S;

[eVecs, Eigenvalues] = eig(A_tilde);

%Gets the DMD eigenvectors back
Eigenvectors=Y*V*inv(S)*eVecs;

Ppf = sum(abs(Eigenvectors),2);