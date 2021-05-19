function [mod,B,V,dev,stats] = alogisticreg(X,g)

if any(g==0)
    g = g + 1;
end

% remove minimum variance params fot stability
good = find(var(X) > 1e-2);

V = zeros(size(X,2),length(good));
for i = 1:length(good)
    V(good(i),i) = 1;
end

% run the logreg
[B,dev,stats] = mnrfit(X(:,good),g);

% reshape the outputs
mod.c = B(1);
mod.b = B(2:end)'*V';

mod.t = stats.t(2:end)'*V';
mod.p = stats.p(2:end)'*V';
mod.covb = V*stats.covb(2:end,2:end)*V';