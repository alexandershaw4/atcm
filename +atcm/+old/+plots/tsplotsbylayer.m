function tsplotsbylayer(s,pst)
% membrane potential timeseries plots by layer

c = { 'k' , 'r' }; % condition colours
n = {'L2/3 (SP)' 'L2/3 (SI)' 'L4 (SS)' 'L5 (DP)' 'L5 (DI)' 'L6 (TP)' 'Thal (RT)' 'Thal (RL)'};

%        sp si ss dp di tp rt rl
order = [2  3  1  4  5  6  7  8]; 

for node = 1:size(s{1},1) % new fig for each 'region' (node)
    figure('position',[1341          40        1016         945]);
    
    
    for i  = 1:length(s) % diff lines for each condition / trial type
        si = squeeze(s{i}(node,:,1,:));

        for j = 1:8
            subplot(8,1,j), plot(pst,si(order(j),:),c{i} ); hold on;
            ylabel('mV');
            xlabel(sprintf('%s - time (ms)',n{j}));
        end
    end

    set(findall(gcf,'-property','FontSize'),'FontSize',14)
end