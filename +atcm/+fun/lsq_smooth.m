function [y,model] = lsq_smooth(x)
% Non-negative least squares smoothing of input vector x.
%
%   [y,model] = atcm.fun.lsq_smooth(x)
%
% Assumes each point in sequence x could be represented by a Gaussian with
% fixed variance, then fits this Gaussian mixture back to the data to find
% which components are needed to explain original sequence.
%
% AS2022

x = x(:);
g = AGenQn(x,6);
b = atcm.fun.lsqnonneg(real(g),real(x));
y = g*b;

% outputs
model.g = g;
model.b = b;
model.x = x;
model.centres = find(b);