function [f,J,Q,D] = tcm_b(x,u,P,M,m)
% State equations for an extended canonical thalamo-cortical neural-mass model.
%
% This model implements a conductance-based canonical thalamo-cortical circuit,
% with cytoarchitecture inspired by Gilbert & Wiesel (1983), Douglas & 
% Martin (2004) and Traub (2004) models.
%
% The equations of motion are Moris Lecar-esque equations, similar to Moran
% (2011), but with conductances for AMPA, NMDA, GABA-A, & GABA-B channels. 
% These 'channels' feature their own reversal poentials and rate constants:
%
% K  = -70           (Leak)
% Na =  60  / 4 ms   (AMPA)
% Cl = -90  / 16 ms  (GABA-A)
% Ca =  10  / 100 ms (NMDA)   + voltage mag switch
% B  = -100 / 200 ms (GABA-B)
% f  = -40
%
% FORMAT [f,J,Q,D] = atcm.tcm_hilge(x,u,P,M)
%
% x - states and covariances
%
% x(i,j,k)        - k-th state of j-th population of i-th source
%                   i.e., running over sources, pop. and states
%
%   population: 1  - Spint stellates (L4)
%               2  - Superficial pyramids (L2/3)
%               3  - Inhibitory interneurons (L2/3)     
%               4  - Deep pyramidal cells (L5)
%               5  - Deep interneurons (L5)
%               6  - Thalamic projection neurons (pyramid) (L6)
%               7  - Reticular cells (Thal)
%               8  - Thalamo-cortical relay cells (Thal)
%
%
%        state: 1 V   - voltage
%               2 gE  - conductance: AMPA   (excitatory)
%               3 gI  -  conductance: GABA-A (inhibitory)
%               4 gN  - conductance: NMDA   (excitatory)
%               5 gB  - conductance: GABA-B (inhibitory)
%               6 gM  - conductance: M-channels (inhibitory)
%               7 gih - conductance: H-channels (inhibitory)
%
%      outputs: f = model states as a vector - hint: spm_unvec(f,M.x) 
%               J = system Jacobian - dfdx
%               Q = delay operator  - Q = inv(1 - D.*dfdx)*f(x(t))
%               D = states delay matrix
%
% Info:
%  - Ih is a hyperpolarization-activated cation current mediated by HCN channel
%  - non-selective, voltag gated, responsible for cariac 'funny' (pacemaker) current
%
%  - M-channels (aka Kv7) are noninactivating potassium channels
%  - M is unique because it is open at rest and even more likely to be open during depolarization
%  - M is a pip2 regulated ion channel
%
% Alexander Shaw 2019: ShawA10@cardiff.ac.uk
%
% Notes, changes, updates:
%
% Extrinsics connection matrices [ampa but AN{n} is nmda equiv]:
% A{1} = Forward  SP -> SS & DP
% A{2} = Backward DP -> SP & SI
% A{3} = Back/Lat TP -> SS & TP
% A{4} = Inter-Thal [B] RT -> RC
% A{5} = Inter-Thal [F] RC -> RT
%--------------------------------------------------------------------------



% Flag: include M- & H- channels on L6 TP & Thalamic Relay cells, or not
%--------------------------------------------------------------------------
IncludeMH = 1;
 
% get dimensions and configure state variables
%--------------------------------------------------------------------------
ns   = size(M.x,1);                      % number of sources
np   = size(M.x,2);                      % number of populations per source
nk   = size(M.x,3);                      % number of states per population
x    = reshape(x,ns,np,nk);              % hidden states 


% extrinsic connection strengths
%==========================================================================
 
% exponential transform to ensure positivity constraints
%--------------------------------------------------------------------------
for i = 1:length( P.A )
    A{i}  = exp(P.A{i});
    AN{i} = exp(P.AN{i});
end

C     = exp(P.C); 
 

% detect and reduce the strength of reciprocal (lateral) connections
%--------------------------------------------------------------------------
for i = 1:length(A)
    L    = (A{i} > exp(-8)) & (A{i}' > exp(-8));
    A{i} = A{i}./(1 + 8*L);
end

            
% intrinsic connection strengths
%==========================================================================
G    = full(P.H);
G    = exp(G);

% nmda matrix
%P.Hn=P.H;
%P.Hn(3,3)=P.nmda;

Gn = full(P.Hn);
Gn = exp(Gn);
%Gn = exp(full(G));

% this was for specifying trial-specific intrinsic connections (w Betas)
% in the LTP project (Sumner, Spriggs, Shaw 2020)
%--------------------------------------------------
% if ~all(size(P.G)==np)
%     for i = 1:size(G,3)
%         % trial specific intrinsic effects !
%         Gtrial   = diag( (P.G));
%         G(:,:,i) = G(:,:,i) + Gtrial; 
%     end
% elseif all(size(P.G)==np)
%     % a full 8*8 connectivity for this trial
%     for i = 1:size(G,3)
%         G(:,:,i) = G(:,:,i) + P.G;
%     end
% end




% connectivity switches
%==========================================================================
% 1 - excitatory spiny stellate cells (granular input cells)
% 2 - superficial pyramidal cells     (forward  output cells)
% 3 - inhibitory interneurons         (intrisic interneuons)
% 4 - deep pyramidal cells            (backward output cells)
% 5 - deep interneurons               
% 6 - thalamic projection pyramidal cells (with m- and h- currents)
% 7 - thalamic reticular cells
% 8 - thalamic relay cells (with m- and h- currents)
%
% Thalamic cells attached to different cortical regions (models) are laterally connected


% % extrinsic connections (F B) - from superficial and deep pyramidal cells
% %--------------------------------------------------------------------------
%       SP  DP  tp  rt  rc
SA   = [1   0   1   0   0;   %  SS    % added TP->SP
        0   1   0   0   0;   %  SP
        0   1   0   0   0;   %  SI
        1   0   0   0   0;   %  DP
        0   0   0   0   0;   %  DI
        0   0   1   0   0;   %  TP % 0 in ket study
        0   0   0   0   1;   %  rt % 0 in ket study
        0   0   0   1   0]/8;%  rc % 0 in ket study
    
    SA(:,[3 4 5]) = 0; % For ket study
    
    
% % extrinsic NMDA-mediated connections (F B) - from superficial and deep pyramidal cells
% %--------------------------------------------------------------------------    
SNMDA = [1   0   1   0   0;   %  SS
         0   1   0   0   0;   %  SP
         0   1   0   0   0;   %  SI
         1   0   0   0   0;   %  DP
         0   0   0   0   0;   %  DI
         0   0   1   0   0;   %  TP % 0 in ket study
         0   0   0   0   1;   %  rt % 0 in ket study
         0   0   0   1   0]/8;%  rc % 0 in ket study

     SNMDA(:,[3 4 5]) = 0; % For ket study
     
% intrinsic connectivity switches
%--------------------------------------------------------------------------    
%   population: 1  - Spint stellates (L4)                : e
%               2  - Superficial pyramids (L2/3)         : e
%               3  - Inhibitory interneurons (L2/3)      : i
%               4  - Deep pyramidal cells (L5)           : e
%               5  - Deep interneurons (L5)              : i
%               6  - Thalamic projection neurons -L6     : e
%               7  - Reticular cells (Thal)              : i
%               8  - Thalamo-cortical relay cells (Thal) : e

GEa = zeros(8,8);
GIa = zeros(8,8);

% Excitatory (np x np): AMPA & NMDA
%--------------------------------------------------------------------------
% This is a simplified, predictive-coding friendly excitatory architecture
% NEW - BASED ON HILGETAG PAPER
%           ss  sp  si  dp  di  tp  rt  rl   
GEa(1,:) = [0   4   0   0   0   2   0   2];
GEa(2,:) = [6   4   0   0   0   0   0   1];
GEa(3,:) = [0   8   0   4   0   4   0   0]; 
GEa(4,:) = [1   6   0   0   0   2   0   1];
GEa(5,:) = [0   4   0   4   0   2   0   0];
GEa(6,:) = [0   0   0   6   0   0   0   1/4]; 
GEa(7,:) = [0   0   0   1   0   2   0   2]; 
GEa(8,:) = [0   0   0   1   0   2   0   0];

% GEa(1,:) = [0   0   0   0   0   0   0   2]/1;
% GEa(2,:) = [4   0   0   0   0   0   0   0]/1;
% GEa(3,:) = [0   4   0   0   0   0   0   0]/1; 
% GEa(4,:) = [0   4   0   0   0   0   0   0]/1;
% GEa(5,:) = [0   0   0   2   0   2   0   0]/1;
% GEa(6,:) = [0   0   0   2   0   0   0   1]/1; 
% GEa(7,:) = [0   0   0   0   0   0   0   2]/1; 
% GEa(8,:) = [0   0   0   0   0   2   0   0]/1;

% added TP>SI = 2 (was 0)
GEa = GEa./8;
%GEa = GEa .* ~eye(np);
GEa = GEa + eye(np);      % KILLED

GEn = GEa;
GEn = GEn + eye(8)/8;

% Inhibitory connections (np x np): GABA-A & GABA-B
%--------------------------------------------------------------------------
%           ss  sp  si  dp  di  tp  rt  rl
GIa(1,:) = [16  0   8   0   0   0   0   0 ];
GIa(2,:) = [0   32  16  0   0   0   0   0 ]; %spsp was 16
GIa(3,:) = [0   0   32  0   0   0   0   0 ];
GIa(4,:) = [0   0   0   8   12  0   0   0 ];
GIa(5,:) = [0   0   32  0   16  0   0   0 ];
GIa(6,:) = [0   0   0   0   32  8   0   0 ];
GIa(7,:) = [0   0   0   0   0   0   32  0 ];
GIa(8,:) = [0   0   0   0   0   0   8   32]; 

% GIa(1,:) = [2   0   16  0   0   0   0   0 ];
% GIa(2,:) = [0   64  16  0   0   0   0   0 ]; %spsp was 16
% GIa(3,:) = [0   0   32  0   0   0   0   0 ];
% GIa(4,:) = [0   0   0  128  12  0   0   0 ];
% GIa(5,:) = [0   0   32  0   16  0   0   0 ];
% GIa(6,:) = [0   0   0   0   32  8   0   0 ];
% GIa(7,:) = [0   0   0   0   0   0   32  0 ];
% GIa(8,:) = [0   0   0   0   0   0   8   32]; 

%GIa(2,3)=54;
%GIa(3,3)=128+64;
%GIa(5,5)=64;
%GIa(6,6)=16;

GIa = GIa/2;
%GIa = GIa/1.2;

%GIa = 6*~~GIa;

GIb = GIa;

if isfield(P,'scale')
    GEa = GEa * exp(P.scale(1));
    GEn = GEn * exp(P.scale(2));
    GIa = GIa * exp(P.scale(3));
    GIb = GIb * exp(P.scale(4));
end


if IncludeMH
    
    % M- & H- channel conductances (np x np) {L6 & Thal Relay cells only}
    %----------------------------------------------------------------------
    VM   = -70;                            % reversal potential m-channels          
    VH   = -30;                            % reversal potential h-channels 

    GIm  = sparse([6 8],[6 8],4,8,8);
    Mh   = diag(exp(P.Mh));

    GIh      = full(sparse([6 8],[6 8],4   ,8,8)); % 1/4
    Hh       = exp(P.Hh);
    GIh(6,6) = GIh(6,6)*Hh(1);
    GIh(8,8) = GIh(8,8)*Hh(2);

    KM    = (exp(-P.m)*1000/160) ;               % m-current opening + CV
    KH    = (exp(-P.h)*1000/100) ;               % h-current opening + CV
    h     = 1 - spm_Ncdf_jdw(x(:,:,1),-100,300); % mean firing for h-currents
end

% Channel rate constants [decay times]
%--------------------------------------------------------------------------
KE  = exp(-P.T(:,1))*1000/4;            % excitatory rate constants (AMPA)
KI  = exp(-P.T(:,2))*1000/16;           % inhibitory rate constants (GABAa)
KN  = exp(-P.T(:,3))*1000/100;          % excitatory rate constants (NMDA)
KB  = exp(-P.T(:,4))*1000/200;          % excitatory rate constants (NMDA)

% now using faster AMPA and GABA-A dynamics based on this book:
% https://neuronaldynamics.epfl.ch/online/Ch3.S1.html#:~:text=GABAA%20synapses%20have%20a,been%20deemed%203%20times%20larger.

%KE  = exp(-P.T(:,1))*1000/3;            % excitatory rate constants (AMPA)
%KN  = exp(-P.T(:,3))*1000/150;          % excitatory rate constants (NMDA)
KI  = exp(-P.T(:,2))*1000/6;           % inhibitory rate constants (GABAa)

% Trial effects on time constants: AMPA & NMDA only
if isfield(P,'T1')
    KE = KE + P.T1(1);
    KN = KN + P.T1(2);
end

% Voltages [reversal potentials] (mV)
%--------------------------------------------------------------------------
VL   = -70;                               % reversal  potential leak (K)
VE   =  60;                               % reversal  potential excite (Na)
VI   = -90;                               % reversal  potential inhib (Cl)
VR   = -40;   %55                            % threshold potential (firing)
VN   =  10;                               % reversal Ca(NMDA)   
VB   = -100;                              % reversal of GABA-B

% membrane capacitances {ss  sp  ii  dp  di  tp   rt  rl}
%--------------------------------------------------------------------------
CV   = exp(P.CV).*      [128 128 128  128 64  128  64  64*2]/1000;  

CV   = exp(P.CV).*      [128 32  32  128 64  128  256 64]/1000;  


% [128 128 256 32]

% leak conductance - fixed
%--------------------------------------------------------------------------
GL   = 1 ;       

% mean-field effects:
%==========================================================================

% neural-mass approximation to covariance of states: trial specific
%----------------------------------------------------------------------
Vx   = exp(P.S)*32; % 32
%    {ss    sp    ii    dp    di    tp    rt    rl}
%Vx = [Vx(1) Vx(1) Vx(2) Vx(1) Vx(2) Vx(1) Vx(3) Vx(3)];

if nargin < 5
    % compute only if not passed by integrator
    m    =     spm_Ncdf_jdw(x(:,:,1),VR,Vx);
    
end

% extrinsic effects
%--------------------------------------------------------------------------
a       = zeros(ns,5);
an      = zeros(ns,5); 
a(:,1)  = A{1}*m(:,2);                      % forward afference  AMPA - SP->SS
a(:,2)  = A{2}*m(:,4);                      % backward afference AMPA - DP->SP
a(:,3)  = A{3}*m(:,6);                      % FWD thalamic projection pyramids
a(:,4)  = A{4}*m(:,7);                      % LAT reticular AMPA
a(:,5)  = A{5}*m(:,8);                      % LAT relay AMPA
an(:,1) = AN{1}*m(:,2);                     % forward afference  NMDA
an(:,2) = AN{2}*m(:,4);                     % backward afference NMDA
an(:,3) = AN{3}*m(:,6);                     % thalamic projection pyramids
an(:,4) = AN{4}*m(:,7);                     % reticular NMDA
an(:,5) = AN{5}*m(:,8);                     % relay NMDA

% Averge background activity and exogenous input
%==========================================================================
BE     = exp(P.E)*0.8;

% input(s)
%--------------------------------------------------------------------------
if isfield(M,'u')
      U =   u;%(:); % endogenous input
else; U = C*u;%(:); % exogenous input
end

% flow over every (ns x np) subpopulation
%==========================================================================
f     = x;

% Thalamo-cortical flow [eq. motion] over modes, populations, states...
%--------------------------------------------------------------------------
for i = 1:ns
   
        % allow switching of thalamus on/off
        %HasThal = M.HasThal(i);
            
        % input scaling: 
        %------------------------------------------------------------------
        if any(full(U(:))) ;
            dU = u(1)*C(i,1) ;
            %dU = u(1)*( C(i,:).*[1 1/64 1/64 1/64 1/64] );
        else
            dU = 0;
        end
        
        Gsc = ~eye(8);
        Gsc = Gsc +  (diag(exp(P.Gsc)));
        
        % intrinsic coupling - parameterised
        %------------------------------------------------------------------
        E      = ( (G(:,:,i).*GEa).*Gsc )*m(i,:)'; % AMPA currents
        ENMDA  = (Gn(:,:,i).*GEn)*m(i,:)'; % NMDA currents
        I      = ( G(:,:,i).*GIa)*m(i,:)'; % GABA-A currents
        IB     = ( G(:,:,i).*GIb)*m(i,:)'; % GABA-B currents
        
        if IncludeMH
            
            % intrinsic coupling - non-parameterised: intrinsic dynamics
            %--------------------------------------------------------------
            Im     = (Mh(:,:).*GIm)*m(i,:)'; % M currents
            Ih     =           GIh *h(i,:)'; % H currents
        end
        
        % extrinsic coupling (excitatory only) and background activity
        %------------------------------------------------------------------
        E     = (E     +  BE  + SA   *a (i,:)')*2;
        ENMDA = (ENMDA +  BE  + SNMDA*an(i,:)')*2;
      
        % and exogenous input(U): 
        %------------------------------------------------------------------
        input_cell        = 8;%[8 1 2 4 6];
        
        if isfield(M,'inputcell');
            input_cell = M.inputcell;
        end
        
        %dU = dU*[1 1 0 0 0];
        
        E(input_cell)     = E(input_cell)         +dU';
        ENMDA(input_cell) = ENMDA(input_cell)     +dU';
                
        % Voltage equation
        %==================================================================
        if ~IncludeMH
            
          f(i,:,1) =         (GL*(VL - x(i,:,1))+...
                       1.0*x(i,:,2).*(VE - x(i,:,1))+...
                       1.0*x(i,:,3).*(VI - x(i,:,1))+...
                       1.0*x(i,:,5).*(VB - x(i,:,1))+...
                       1.0*x(i,:,4).*(VN - x(i,:,1)).*mg_switch(x(i,:,1)))./CV;
            
        elseif IncludeMH
            
          f(i,:,1) =         (GL*(VL - x(i,:,1))+...
                       x(i,:,2).*(VE - x(i,:,1))+...
                       x(i,:,3).*(VI - x(i,:,1))+...
                       x(i,:,5).*(VB - x(i,:,1))+...
                       x(i,:,6).*(VM - x(i,:,1))+...
                       x(i,:,7).*(VH - x(i,:,1))+...
                       x(i,:,4).*(VN - x(i,:,1)).*mg_switch(x(i,:,1)))./CV;
        end
                   
        % Conductance equations
        %==================================================================   
        %pop_rates = [1 1 2 1 1 1 1/8 1/8];
        pop_rates = [1 1 1 1 1 1 1 1];
        %pop_rates = pop_rates.*exp(P.pr);
        %gabaa_rate = pop_rates .* exp(P.gaba);
        
        %if isfield(P,'nmdat')
        %    nmdat = P.nmdat;
        %else
        %    nmdat = pop_rates;
        %end
        
        f(i,:,2) = (E'     - x(i,:,2)).* (KE(:,:)');%*pop_rates);
        f(i,:,3) = (I'     - x(i,:,3)).* (KI(:,:)');%*gabaa_rate);
        f(i,:,5) = (IB'    - x(i,:,5)).* (KB(:,:)');%*pop_rates);
        f(i,:,4) = (ENMDA' - x(i,:,4)).* (KN(:,:)');%*nmdat);
        
        if IncludeMH
            f(i,:,6) = (Im'    - x(i,:,6)).*(KM(i,:) );%*pop_rates );
            f(i,:,7) = (Ih'    - x(i,:,7)).*(KH(i,:) );%*pop_rates );
        end
        
        % c.f. synaptic delays + conduction delays
        %------------------------------------------------------------------
        DV = 1./[1 1 1 1 1 1 1 1]; % this is usually switched off
        
        DV       = 1./[1 1 1 2.2 1 2 8 8]; % intrinsc are already 10ms, so 10*8

        
        if isfield(P,'TV')
            DV       = DV.*exp(P.TV);
            f(i,:,2) = f(i,:,2) .* DV;  % AMPA
            f(i,:,3) = f(i,:,3) .* DV;  % GABA-A
            f(i,:,4) = f(i,:,4) .* DV;  % NMDA
            f(i,:,5) = f(i,:,5) .* DV;  % GABA-B

            if IncludeMH
                f(i,:,6) = f(i,:,6) .* DV;  % M
                f(i,:,7) = f(i,:,7) .* DV;  % H
            end 
        end
        
                
end


% vectorise equations of motion
%==========================================================================
f = spm_vec(f);
 
if nargout < 2, return, end

% Only compute Jacobian (gradients) if requested
%==========================================================================
J = spm_cat(spm_diff(M.f,x,u,P,M,1));

if nargout < 3, return, end

% Only compute Delays if requested
%==========================================================================
% Delay differential equations can be integrated efficiently (but 
% approximately) by absorbing the delay operator into the Jacobian
%
%    dx(t)/dt     = f(x(t - d))
%                 = Q(d)f(x(t))
%
%    J(d)         = Q(d)df/dx
%--------------------------------------------------------------------------
% [specified] fixed parameters
%--------------------------------------------------------------------------
D  = [1 16];
d  = -D.*full(exp(P.D(1:2)))/1000;
Sp = kron(ones(nk,nk),kron( eye(np,np),eye(ns,ns)));  % states: same pop.
Ss = kron(ones(nk,nk),kron(ones(np,np),eye(ns,ns)));  % states: same source

% Thalamo cortical interactions: ~80ms round trip: 20 ms T->C, 60 ms C->T
%--------------------------------------------------------------------------
%Thalamocortical connections and forward connections from Vp to Vs had
%a mean delay of 3 ms, while corticothalamic connections and backward
%connections from Vs to Vp had a mean delay of 8 m - Neural Dynamics in a Model of the
%Thalamocortical System. I. Layers, Loops and the Emergence of Fast Synchronous Rhythms
% Lumer et al 1997

CT = 8; %60;
TC = 3; %20;

%CT = 60;
%TC = 20;

Tc              = zeros(np,np);
Tc([7 8],[1:6]) = CT  * exp(P.D0(1)); % L6->thal
Tc([1:6],[7 8]) = TC  * exp(P.D0(2)); % thal->ss
%Tc = Tc.*~~(GEa | GIa);
Tc = -Tc / 1000;
Tc = kron(ones(nk,nk),kron(Tc,eye(ns,ns)));

% Mean intra-population delays, inc. axonal etc. Seem to help oscillation
%--------------------------------------------------------------------------
Dp = ~Ss;                            % states: different sources
Ds = ~Sp & Ss;                       % states: same source different pop.
%Ds = Ds.*(~(Ds & Tc));              % remove t-c and c-t from intrinsic

if ~isfield(P,'delays')
    D  = d(2)*Dp + d(1)*Ds + Tc  ;
else
    D = d(1)*Ds + Tc  ;       %+ Dself;% Complete delay matrix
end

%D = d(2)*Dp + Tc; %%%%%!!!!!!

% Implement: dx(t)/dt = f(x(t - d)) = inv(1 - D.*dfdx)*f(x(t))
%                     = Q*f = Q*J*x(t)
%--------------------------------------------------------------------------
Q  = spm_inv(speye(length(J)) - D.*J);
