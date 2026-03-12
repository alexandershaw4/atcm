function G = spectralEFE(a, Phi, DCM)
    % Get spectral input
    uomega = Phi * a(:);  % action: shape of spectral input
    DCM.M.external_spectrum = uomega;

    % Forward model
    [y,~,~,~] = atcm.fun.Alex_LaplaceTFwD(spm_unvec(DCM.Ep, DCM.M.pE), DCM.M, DCM.xU);

    % Observation prediction error
    err = DCM.xY.y{:} - y{:};
    err = err(:);

    G = sum(err.^2);  % could add entropy/epistemic terms
end