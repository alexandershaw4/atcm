

gV    = 1:8; 
gAMPA = 9:16;
gNMDA = 25:32;

freqs = 1:120;
[Sy, Sx] = atcm.fun.tcm_psd_from_jacobian(J, 1, eye(size(J,1)), 0, freqs); % state PSDs
% Average group power & cross-power
Saa = squeeze(mean(mean(Sx(gAMPA, gAMPA, :),1),2));
Snn = squeeze(mean(mean(Sx(gNMDA, gNMDA, :),1),2));
San = squeeze(mean(mean(Sx(gAMPA, gNMDA, :),1),2));
Coh_AN = abs(San).^2 ./ (Saa .* Snn + eps);

figure; plot(freqs, Coh_AN); xlabel('Hz'); ylabel('Coherence (AMPA–NMDA)');
title('AMPA–NMDA coherence from linearisation');