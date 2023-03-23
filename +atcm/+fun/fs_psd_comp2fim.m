function FIM = fs_psd_comp2fim(x,w,N)
% Compute the approximate (N-dim) fisher im for a psd vector 
% by decomposition to a Gaussianmixture, PCA into components and taking the
% inner product. Taking the trace of this will further estimate the 
% Hutchinson's estimator
%
% AS23

[V] = atcm.fun.GaussAR(x,w,1,N,1);

FIM = V'*V;