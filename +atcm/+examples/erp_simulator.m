% use the thalamocortical model as a MEG single channel ERP simulator


% function handles
M.IS = @atcm.integrate_erp;
M.f  = @atcm.tc_dev_dev;

% time, freqs etc.
M.Hz      = 4:80;
dt        = 1/1200;
M.sim.pst = (0:dt:1)*1000;
M.sim.dt  = dt;
    
% initial states
M.x = zeros(1,8,7);
M.x(:,:,1) = -70;

DCM.M = M;

% neural priors
DCM = atcm.parameters(DCM,1);

% simulus / input bump
R(1) = 0.69; % input bump (stim) delay: exp(R(1))*60
R(2) = .5;   % input bump (stim) size:  exp(R(2))*8
R(3) = .1;    % input bump (stim) width
DCM.M.pE.R = R;

% states that contirbute to the signal
J = zeros(1,56)-1000;   % 56 total per node
J([1 2 4 6]) = log(.8); % pyramidal cells only
DCM.M.pE.J = J;

% Trial data
U.X = [0];

% generate data
for j = 1:2 % standard / deviant
    for i = 1:40

        R0(1) = j*2*R(1) + randn(1);
        R0(2) = j*1*R(2) + randn(1); % in condition 1 the input is on average 2*bigger
        R0(3) = 2 + ( .2*randn(1) );
        
        % Generate the stimulus time-course
        d = 60 * R0(1);
        s = 8  * exp(R0(2));
        w = 1 * exp(R0(3));
        drive(j,i,:)  = atcm.fun.makef(M.sim.pst,d,s,w);
        
        P     = DCM.M.pE;
        P.R   = R0;

        [y,w,s,g,t,pst] = feval(DCM.M.IS,P,DCM.M,U);

        data{j}(i,:) = g{1};
        spec{j}(i,:) = s{1};
    end
end

% plot the model simulations for each trial type
%--------------------------------------------------------------------------
%scale(1) = ( min([ min(data{1}(:)) min(data{2}(:)) ]) );
%scale(2) = ( max([ max(data{1}(:)) max(data{2}(:)) ]) );

figure,

subplot(3,2,[1 2]);
plot(pst,squeeze(mean(drive,2)));
legend({'Stim 1' 'Stim 2'});

ds = 200;
dpst = pst(1:ds:end);
dpst = round(dpst/10)*10;
subplot(3,2,3);imagesc(real(data{1})); %caxis([scale]);
ylabel('trial'); xlabel('time');
set(gca,'xtick',1:ds:length(pst),'xticklabels',dpst);
title('condition 1');

subplot(3,2,4);imagesc(real(data{2})); %caxis([scale]);
ylabel('trial'); xlabel('time');
set(gca,'xtick',1:ds:length(pst),'xticklabels',dpst);
title('condition 2');

subplot(3,2,[5 6]);
plot(pst,mean(data{1},1),pst,mean(data{2},1));
legend({'Mean condition 1' , 'Mean condition 2'});
xlabel('time'); ylabel('voltage (mV)');
set(findall(gcf,'-property','FontSize'),'FontSize',20)