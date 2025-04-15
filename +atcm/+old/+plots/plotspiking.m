function plotspiking(sp,pst)
% plots firing over time!
%
% input the cell array 'Sp' returned by the integrator, e.g.
%
% [y,w,s,g,t,pst,l,n,f,Qd,Sp] = feval(DCM.M.IS,DCM.Ep,DCM.M,DCM.xU);
%
% atcm.plots.plotspiking(Sp,pst)
%
%

if iscell(sp)
    ntr = length(sp);
else
    ntr = 1;
    sp = {sp};
end

for nt = 1:length(ntr)
    figure; % new figure for each trial type in model
    
    % plot 1: cells firing over time
    f = squeeze(sp{nt});
    subplot(121); plot(pst(f(:,1)),f(:,2),'.','MarkerSize',15)
    cells = {'ss' 'sp' 'si' 'dp' 'di' 'tp' 'rt' 'rl'};
    set(gca,'ytick',1:8,'yticklabels',cells);
    xlabel('time (ms)');
    title('firing');
    
    % plot 2: histogram
    subplot(122);
    hist3(f,[10 8]); % cells by time, chunked into 10 time bins
    set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
    
    
end

    
set(findall(gcf,'-property','FontSize'),'FontSize',20)