
if ~exist('DCM','var') && exist('TCM','var')
    DCM = TCM;
    flagremoveDCM = 1;
else
    flagremoveDCM = 0;
end

% generate tcm diagnostic plots
%figure('position',[1010          42         774         943]);
figure

r2 = 100*corr( real(DCM.xY.y{1}(:)), real(y{1}(:)) ).^2;
r2 = round(r2);
subplot(221); plot(w,DCM.xY.y{1},':',w,y{1});
title(sprintf('Data & Model - Fit\n(r^2 = %d%%)',r2));
xlabel('Frequency (Hz)');ylabel('PSD');
grid on;

try
subplot(222); plot(w,squeeze(l{1}.iweighted(1,:,:)),'linewidth',2);
title('Principal (contributing) cells spectra');
xlabel('Frequency (Hz)');ylabel('PSD');
grid on;
end

subplot(2,2,[3 4]); plot(pst,squeeze(s{1}(1,:,1,:)),'linewidth',2);
title('Membrane potential time series');
xlabel('Time (s)');ylabel('mV');
grid on;

set(findall(gcf,'-property','FontSize'),'FontSize',20)

if flagremoveDCM
    clear DCM flagremoveDCM;
end