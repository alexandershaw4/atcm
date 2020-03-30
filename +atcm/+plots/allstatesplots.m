function allstatesplots(s,pst)

ntr = length(s);
[ns,npp,nk,nt] = size(s{1});

staten = {'mV' 'AMPA' 'GABA A' 'NMDA' 'GABA B' 'M' 'H'};
pops   = {'ss' 'sp' 'si' 'dp' 'di' 'tp' 'rt' 'rl'};
    
for i = 1:ntr
    for is = 1:ns
        for j = 1:nk
            figure;
            data = squeeze(s{i}(is,:,j,:));
            for p = 1:npp
                subplot(2,4,p); plot(pst,data(p,:),'b','linewidth',2);
                title(sprintf('%s [%s]',upper(pops{p}),upper(staten{j})));
            end
            drawnow;
        end
    end
end