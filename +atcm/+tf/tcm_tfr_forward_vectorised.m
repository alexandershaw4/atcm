function [y_vec, cache] = tcm_tfr_forward_vectorised(P, M, U, tfopts, data)
% Wraps forward model and returns a vector for VL
% data.TFR: nf × nt (target % change), data.W: nf × nt weights (0..1)
[TFRm, t, f] = atcm.tf.tcm_timefreq_tf(P,M,U,tfopts);

% Optional transforms to stabilise the distribution
X  = asinh(TFRm/50);              % variance stabilisation (robust)
Xd = asinh(data.TFR/50);

% Weights (e.g., emphasise 30–90 Hz and 0.1–1.0 s)
if ~isfield(data,'W') || isempty(data.W)
    W = ones(size(X));
else
    W = data.W;
end

% Compute residual matrix (for diagnostics)
R = W .* (X - Xd);

% Vectorise for your fitter
y_vec = R(:);

% Return cache if you want gradients / plotting
cache = struct('t',t,'f',f,'TFRm',TFRm,'X',X,'Xd',Xd,'W',W);
end
