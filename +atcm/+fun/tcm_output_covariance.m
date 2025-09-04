function Sigma_y = tcm_output_covariance(Sigma_x, C, R)
% Maps state covariance to observation space: Σ_y = C Σ_x C' + R
    if isscalar(R), R = R * eye(size(C,1)); end
    Sigma_y = C * Sigma_x * C' + R;
end
