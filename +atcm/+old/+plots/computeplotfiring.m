function computeplotfiring(f,dt,pst)
% pass f, dt & pst. f (& pst) returned by integration scheme atcm.integrate
%
% computes firing rate over time per population - for each region and trial
% type in the model.
%
% 

c = { 'k' , 'r' }; % condition colours
n = {'L2/3 (SP)' 'L2/3 (SI)' 'L4 (SS)' 'L5 (DP)' 'L5 (DI)' 'L6 (TP)' 'Thal (RT)' 'Thal (RL)'};

%        sp si ss dp di tp rt rl
order = [2  3  1  4  5  6  7  8]; 

% Compute firing rate from f returned by integration scheme, using f2sr.m
for i = 1:length(f) % per trial
    f0(i,:,:,:) = atcm.fun.f2sr(f{i},dt);
end

for i = 1:size(f0,2) % regions modelled
    figure('position',[1224          40         588         945]);
    for j = 1:size(f0,1) % trials
        for k = 1:8
            subplot(8,1,k), plot(pst, squeeze(f0(i,j,order(k),:)),c{i} ); hold on;
            ylabel('Firing rate (Hz)');
            xlabel(sprintf('%s - time (ms)',n{k}));
        end
    end
end

