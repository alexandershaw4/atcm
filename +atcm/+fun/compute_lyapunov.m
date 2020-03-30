function [LLE, lambda] = compute_lyapunov(ts,tau,m)
% wrapper on lyaprosen.m to compute Lyapunov exponents using Rosenstein's
% method
%
% ts is chan*time
% tau = delay
% m = dimension [3]
%
% AS

warning off; 

for i = 1:size(ts,1)
    data = ts(i,:);
    try
        [LLE(i) lambda(:,i)] = atcm.fun.lyaprosen(data,tau,m);
    catch
        LLE(i) = 0;
        lambda(:,i) = 0;
        fprintf('LYAP failed for this dataset\n');
    end
end

warning on;
