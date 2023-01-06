function Hz = firing_tcm(P,M,U)
% Compute firing rate (in Hz) from tcm model being integrated with atcm
% integration functions (e.g. atcm.integrate1.m)
%
% Hz = atcm.fun.firing_tcm(P,M,U)
%
% AS2023

[y,w,s,g,drive,pst,l,oth] = feval(M.IS,P,M,U);

% extract spike record from integrator
F = oth.Spike;
t = F{1}(:,1);
f = F{1}(:,2);
S = zeros(8,t(end));

% turn into matrix timseries representation
for i = 1:length(t)
    x = find(t == i);
    n = f(x);
    
    if any(n)
        S(n,i) = 1;
    end
    
end

% divide by length of integrated eriod
tt = [0; diff(pst)]+pst;
tt = tt./1000;
Hz = sum(S,2)./tt(end);