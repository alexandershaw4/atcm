function [y,w] = KuramotoTCM(P,M,varargin)
% these are the state equations as well as the integration and transfer
% functions for the kuramoto thalamo-cortical model
%
%
% P.H = zeros(8,8);
% P.D = zeros(8,8);
% P.fq = zeros(1,8);
% P.k = 0;
% P.d = 0;
% P.J = zeros(1,8)-1000;
% P.J([1 2 4 6]) = log([.2 .8 .2 .2]);
% P.L = 0;
% 
% M.dt    = 1/1200;                   
% M.Fs    = 1/M.dt;                     
% M.tn    = 2;                            
% M.pst   = 1000*((0:M.dt:M.tn-M.dt)');
% M.w     = 4:80; 
% 




C = [...
...  ss    sp    si    dp    di    tp    rt    rl    
     1     0     1     0     0     1     0     1;
     1     1     1     0     0     0     0     0;
     1     1     1     0     0     0     0     0;
     0     1     0     1     1     0     0     0;
     0     0     0     1     1     0     0     0;
     0     0     0     1     1     1     0     1/4;
     0     0     0     0     0     0     1     1;
     0     0     0     0     0     1     1     1 ]/100;


D = [...
   1   0   1   0   0   1   0   300;
   1   1   1   0   0   0   0   0;
   1   1   1   0   0   0   0   0;
   0   1   0   1   1   0   0   0;
   0   0   0   1   1   0   0   0;
   0   0   0   1   1   1   0   300;
   0   0   0   0   0   0   1   1;
   0   0   0   0   0   800   1   1 ]/1000;

C = C.*exp(P.H); % connection strengths
D = C.*exp(P.D); % connection delays

% dist of frequencies:       ss sp si dp di tp rt rl
f_dist                  = ( [10 60 80 20 40 20 10 10].*exp(P.fq) )';


% Parameters 
%------------------------------------------------------------------------
frequency_mean = 10; % mean in Hz per node
f_std          = 0;          % std of natural freqs
k              = exp(P.k);    % global coupling strength
tau            = exp(P.d)*.01;   % mean delay between brain areas
t_max          = M.tn;         % time (s)
dt             = M.dt;    % time-step for integration (1/fs)
sampling       = 1;        % downsampling
sig_n          = 0;        % stdev of noise


% Cabral code
%------------------------------------------------------------------------
[Phases] = atcm.fun.Network_Kuramoto(C,D,frequency_mean,f_std,f_dist,k,tau,t_max,dt,sampling, sig_n);
%Waves    = sin(Phases);
Waves = Phases;

% LFP is a linear mixture of principal cells, the amplitude (mix) of which
% is controlled by J
for i = 1:8
    LFP(i,:) = exp(P.J(i))*Waves(i,:);
end

% % Dynamic mode decomposition to identify the frequency basis set - i.e. the
% % principal frequencies among the contributing populations
% NDMD = 4;
% [Eigenvalues, Eigenvectors, ModeAmplitudes, ModeFrequencies, ...
%         GrowthRates, POD_Mode_Energies] = atcm.fun.dmd(LFP, NDMD, dt);

Ji = find( exp(P.J) );
for i = 1:length(Ji);%NDMD
    
    %[Pf(i,:),Hz]  = atcm.fun.AfftSmooth(Eigenvectors(i,:),1/dt,M.w,60);
    [Pf(i,:),Hz]  = atcm.fun.AfftSmooth(LFP(Ji(i),:),1/dt,M.w,60);
    
end

y{1} = exp(P.L)*double( sum(Pf,1) );
w = M.w;





