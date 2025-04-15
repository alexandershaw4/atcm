function adjplot(h)

p = imagesc(h); axis square;

s = max(abs(h(:)));
caxis([-s s]);
colormap(cmocean('balance'));
colorbar

pop = {'ss' 'sp' 'si' 'dp' 'di' 'tp' 'rt' 'rl'};
set(gca,'xtick',1:8,'xticklabels',pop);
set(gca,'ytick',1:8,'yticklabels',pop);
set(findall(gcf,'-property','FontSize'),'FontSize',16)