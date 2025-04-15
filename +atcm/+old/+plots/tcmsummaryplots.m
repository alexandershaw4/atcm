function tcmsummaryplots(DCM,Ep)

if isvector(Ep)
    Ep = spm_unvec(Ep,DCM.M.pE);
end

[y,w,s,g,drive,pst,l,oth] = feval(DCM.M.IS,Ep,DCM.M,DCM.xU);


% big fig!
figure('position',[515         148        1485        1189]);

% spectra
subplot(8,3,[1 2 4 5 7 8 10 11]);
r2 = 100*corr( real(DCM.xY.y{1}), real(y{1}) ).^2;
r2 = round(r2);
plot(w,DCM.xY.y{1},':',w,y{1});
title(sprintf('Data & Model - Fit\n(r^2 = %d%%)',r2));
xlabel('Frequency (Hz)');ylabel('PSD');
grid on;

s0 = squeeze(s{1}(1,:,:,:));
for i = 1:8
    subplot(8,3,3*i);
    plot(pst,squeeze(s0(i,1,:)));
end

subplot(8,3,[13 14 16 17 19 20 22 23])
plot(w,squeeze(l{1}.iweighted));

end
