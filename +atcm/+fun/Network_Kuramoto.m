function [Phases,dts] = Network_Kuramoto(C,D,frequency_mean,f_std,f_dist,k,tau,t_max,dt,sampling, sig_n)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    General code to compute whole brain network simulations 
%    using the Kuramoto model of coupled oscillators with time delays.
%    To search the optimal working point fitting neuroimaging data, 
%    a loop over parameters (mainly k and tau) may be necessary.
%
%    Code used in 
%    Cabral et al.(2011) NeuroImage 57  130-139
%    Cabral et al.(2012) Pharmacopsychiatry  45(S01), S57-S64
%    Cabral et al.(2014) NeuroImage 90, 423-435
% 
%    Joana Cabral, 2010 joana.cabral@psych.ox.ac.uk
%
%
%   Revised July 2018
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Inputs:
%
% C                    = Connectivity matrix of size NxN, where each entry 
%                        C(n,p) corresponds to coupling weight from area n
%                        to area p.
% D                    = Distance matrix of size NxN, where each entry 
%                        D(n,p) corresponds to the connection length between 
%                        n and p, given in milimeters. 
% frequency_mean       = Natural oscillation frequency (in Hz) of each node
%                        i.e. 40Hz or 60Hz (in the gamma range).
% f_std                = Standard deviation of natural frequencies across 
%                        regions. If 0, we assume identical  oscillators.
% f_dist               = Vector to specify the distribution of intrinsic 
%                        frequencies (Nx1). If all equal, f_dist=ones(N,1)
% k                    = Global coupling stregnth scaling all connections
% tau                  = Mean delay between brain areas (all values in the
%                        delay matrix (stepsDelay) are scaled accordingly
% t_max                = Total time of simulated activity (seconds)
% dt                   = Integration step (must be smaller than the delays) 
%                        i.e. 1e-4, in seconds.
% sampling             = Simulated activity is downsampled before saving
%                        (i.e. if sampling = 10 => 10*dt = 1ms
% sig_n                = Standard deviation of noise (can be zero)
% file_save            = Name of the file to save the simulations
%
% Outputs:  
%
% Phases               = Simulated phases of N oscillators over a total 
%                        time defined by t_max at a resolution dt*sampling.
%                        The phases are in radians between 0 and 2*pi
%                       
% Note: This code saves the simulated data and parameters in the end.
% Beware of disc space for saving. Integration cycles may be removed if
% enough RAM is available.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_t_total   = ceil(t_max/dt);       % total number of time steps 300*0.0001 = 3000000
n_ts  = ceil(n_t_total/sampling);   % same with downsampling ex. 300000

% Integration cycles were implemented due to limited working memory space
n_tcycle  = 10000;                     % number of time steps in each cycle
n_tcs = n_tcycle/sampling;             % same with downsampling
n_cycles   = ceil(n_t_total/n_tcycle); % number of cycles ex. 300

N          = size(C,1);
C          = k*(C/mean(C(:)));     % Normalize to compare with theory
C          = dt*C;                 % Scale the coupling strengths with time step
I          = C>0;                  %
d_m        = mean(D(I));           % Mean distance
v          = d_m/tau;              % Conduction velocity (m/s)
stepsDelay = round(D/(v*dt*1e3));  % Matrix of time steps per delay
sig_noise  = sig_n*sqrt(dt);       % Scale noise per time step

f_diff     = f_dist*f_std;           % define intrinsinc node frequencies.
omega_step = 2*pi*frequency_mean*dt; % 0.0251radians if f=40Hz and dt = 0.0001.
omega_diff = 2*pi*f_diff*dt;
omegas     = omega_step+omega_diff;  % Phase increment per time step (Nx1)


n_td = fix(max(stepsDelay(:)))+10; % number of time steps for maximal delays
n_tp = n_td+n_tcycle;              % Time steps in one cycle including time for delays


th   = zeros(N,n_tp,'single');     % initialize phase timeseries for one cycle
Phases  = zeros(N,n_ts,'single');     % initialize phase timeseries to save


% Initialization

    th(:,1) = 2*pi*rand(N,1);
    for n=1:N
        th(n,1:n_td) = th(n,1)+(0:omegas(n):(n_td-1)*omegas(n));
        th(n,1:n_td) = mod(th(n,1:n_td),2*pi);
    end

% Equations integration

%disp(['Running Kuramoto simulations for ' num2str(N) ' coupled oscillators'])
%disp(['Coupling strength = ' num2str(k) ])
%disp(['Mean delay = ' num2str(tau) 'ms' ])

tic

for c = 1:n_cycles
    %disp(['Cycle = ' num2str(c) ' of ' num2str(n_cycles)])
    th(:,n_td+1:n_tp) = 0;

    if c < n_cycles
        n_tpc = n_tp;
    else
        n_tpc = n_t_total-(n_cycles-1)*n_tcycle+n_td; % nr of steps to complete total time
    end
    
    for t = n_td:n_tpc-1
        dth = omegas + sig_noise*randn(N,1);
        for n = 1:N
            for p = 1:N
                if C(n,p)>0
                 dth(n) = dth(n) + C(n,p)*sin(th(p,t-stepsDelay(n,p))-th(n,t));
                end
            end
        end
        th(:,t+1) = mod(th(:,t)+dth,2*pi);
    end
    ni = (c-1)*n_tcs;
    ns = ceil((n_tpc-n_td)/sampling);
    Phases(:,ni+1:ni+ns) = th(:,n_td+1:sampling:n_tpc);
    
    th(:,1:n_td)      = th(:,n_tp-n_td+1:n_tp);
    T_e = toc;
    %disp(['   (elapsed time = ',num2str(T_e),' s)'])
end

%Phases(:,1:20/(dt*sampling))=[];  % remove initial 20 seconds of simulations
dts=dt*sampling;

%save(file_save,'Phases','C','D','frequency_mean','f_std','f_dist','k','tau','dts','sig_n')