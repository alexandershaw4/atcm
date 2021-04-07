function R = wcor(Y,w)

[T, N] = size(Y);                                                             % T: number of observations; N: number of variables
temp = Y - repmat(w' * Y, T, 1);                                              % Remove mean (which is, also, weighted)
temp = temp' * (temp .* repmat(w, 1, N));                                     % Covariance Matrix (which is weighted)
temp = 0.5 * (temp + temp');                                                  % Must be exactly symmetric
R = diag(temp);                                                               % Variances
R = temp ./ sqrt(R * R');                                                     % Matrix of Weighted Correlation Coefficients
