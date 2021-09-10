
if ~exist('DCM','var') && exist('TCM','var')
    DCM = TCM;
    flagremoveDCM = 1;
else
    flagremoveDCM = 0;
end

% generate tcm diagnostic plots
%figure('position',[1010          42         774         943]);
figure


subplot(221); plot(w,DCM.xY.y{1},':',w,y{1});
title('Data & Model - Fit');
xlabel('Frequency (Hz)');ylabel('PSD');
grid on;

subplot(222); plot(w,squeeze(l{1}.weighted(1,:,:)),'linewidth',2);
title('Principal (contributing) cells spectra');
xlabel('Frequency (Hz)');ylabel('PSD');
grid on;

subplot(2,2,[3 4]); plot(pst,squeeze(s{1}(1,:,1,:)),'linewidth',2);
title('Membrane potential time series');
xlabel('Time (s)');ylabel('mV');
grid on;

set(findall(gcf,'-property','FontSize'),'FontSize',20)

if flagremoveDCM
    clear DCM;
end