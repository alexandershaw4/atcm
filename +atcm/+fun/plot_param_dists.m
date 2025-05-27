function plot_param_dists(mu,covs,ip)

d = length(mu{1});

for i = 1:length(mu)

    for j = 1:length(ip)

        m(i,j) = mu{i}(j);
        v(i,j) = covs{i}(j,j);

    end

end

st = sqrt(v);

% compute upper/lower edges of grid
b = [min(m)*.8; max(m)*1.2];

for j = 1:length(ip)
    grid(:,j) = linspace(b(1,j),b(2,j),100);
end

for j = 1:length(ip)
    for i = 1:length(mu)
        g(j,i,:) = atcm.fun.makef(grid(:,j),m(i,j),1,st(i,j));
    end
end

nx = ceil(sqrt(length(ip)));
ny = ceil(length(ip)/nx);

n = 0;
for i = 1:length(ip)
    subplot(ny,nx,i);
    surf(grid(:,i),1:size(m,1),squeeze(g(i,:,:)));
    shading interp;
end