function timeseries_stacked(pst,s);

cells = {'ss'  'sp'  'si'  'dp'  'di'  'tp'  'rt'  'rl'};
yx = squeeze(s{1}(1,:,1,:));

orderd = fliplr([2 3 1 5 4 6 7 8]);

for i = orderd
    yx(i,:) = yx(i,:) - mean(yx(i,:));
    yx(i,:) = yx(i,:) ./ std(yx(i,:));
    %yx(i,:) = yx(i,:) + (i*2);
end

pos = real(yx(:,1));
offset = pos - (2:2:16)';

yx = yx - offset;

plot(pst,yx);

pos = real(yx(:,1));
set(gca,'ytick',pos,'yticklabels',cells(orderd));

set(findall(gcf,'-property','FontSize'),'FontSize',14);