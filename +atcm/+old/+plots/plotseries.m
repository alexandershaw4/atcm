function plotseries(s,pst,node)

if nargin < 3 || isempty(node)
    node = 1;
end

figure('position',[1000          44         847         934]);
s0 = squeeze(s{1}(node,:,:,:));

subplot(311); 
plot(pst,squeeze(s0(1,1,:)),'k','linewidth',2); hold on;
plot(pst,squeeze(s0(2,1,:)),'b','linewidth',2);
plot(pst,squeeze(s0(3,1,:)),'color',[.8 .3 0],'linewidth',2);
legend({'SS' 'SP' 'SI'});hold off;

subplot(312); 
plot(pst,squeeze(s0(4,1,:)),'b','linewidth',2);hold on;
plot(pst,squeeze(s0(5,5,:)),'color',[.8 .3 0],'linewidth',2);
plot(pst,squeeze(s0(6,1,:)),'k','linewidth',2);
legend({'DP' 'DI' 'TP'});hold off;

subplot(313); 
plot(pst,squeeze(s0(7,1,:)),'color',[.8 .3 0],'linewidth',2);hold on;
plot(pst,squeeze(s0(8,1,:)),'k','linewidth',2);
legend({'RT' 'RL'});hold off;

set(findall(gcf,'-property','FontSize'),'FontSize',18)
