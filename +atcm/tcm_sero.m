function [f,J,D] = tcm_sero(x,u,P,M)
% State equations for an extended canonical thalamo-cortical neural-mass model.
%
% This model implements a conductance-based canonical thalamo-cortical circuit,
% with cytoarchitecture inspired by Gilbert & Wiesel (1983), Douglas & 
% Martin (2004) and Traub (2004) models.
%
% The equations of motion are Moris Lecar-esque equations, similar to Moran
% (2011), but with conductances for AMPA, NMDA, GABA-A, & GABA-B, Kv7 (M), 
% HCN (h) and 5HT-2A channels. 
%
% K   = -70           (Leak)
% Na  =  60  | 2.2 ms (AMPA)
% Cl  = -90  | 5 ms   (GABA-A)
% Ca  =  10  | 100 ms (NMDA)   + voltage mag switch
% B   = -100 | 300 ms (GABA-B)
% H   = -30  | 100 ms
% M   = -52  | 160 ms
% 5HT = -10  | 200
% f   = -52
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
%               3 gI  - conductance: GABA-A (inhibitory)
%               4 gN  - conductance: NMDA   (excitatory)
%               5 gB  - conductance: GABA-B (inhibitory)
%               6 gM  - conductance: M-channels (inhibitory)
%               7 gih - conductance: H-channels (inhibitory)
%               8 g5HT2A - conductance of serontonin 5HT2a on L5 PCs
%
%      outputs: f = model states as a vector - hint: spm_unvec(f,M.x) 
%               J = system Jacobian - dfdx
%               D = state-by-state delay matrix
%
% Info:
%  - Ih is a hyperpolarization-activated cation current mediated by HCN channel
%  - non-selective, voltag gated, responsible for cariac 'funny' (pacemaker) current
%  - HCN = Hyperpolarization-Activated Cyclic Nucleotide-Gated Channels
%
%  - M-channels (aka Kv7) are noninactivating potassium channels
%  - M is unique because it is open at rest and even more likely to be open during depolarization
%  - M is a pip2 regulated ion channel
%
% Notes, changes, updates:
%
% Kv7 channels are actually everwhere - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7530275/#ref90
%
% Extrinsics connection matrices [ampa but AN{n} is nmda equiv]:
% A{1} = Forward  SP -> SS & DP
% A{2} = Backward DP -> SP & SI
% A{3} = Back/Lat TP -> SS & TP
% A{4} = Inter-Thal [B] RT -> RC
% A{5} = Inter-Thal [F] RC -> RT
%
% Dr Alexander Shaw | 2020 | alexandershaw4[@]gmail.com

if isstruct(P) && isfield(P,'p')
    P = P.p;
end
 
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

Gn = full(P.Hn);
Gn = exp(Gn);

if isfield(P,'Gb')
    Gb = exp(full(P.Gb));
else
    Gb = G;
end

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
SA   = [1   0   0   0   0;   %  SS    % added TP->SP
        0   1   0   0   0;   %  SP
        0   1   0   0   0;   %  SI
        0   0   0   0   0;   %  DP
        0   0   0   0   0;   %  DI
        0   0   0   0   0;   %  TP % 0 in ket study
        0   0   0   0   1;   %  rt % 0 in ket study
        0   0   0   1   0]/8;%  rc % 0 in ket study
    
    SA(:,[3 4 5]) = 0; % For ket study
    
% % extrinsic NMDA-mediated connections (F B) - from superficial and deep pyramidal cells
% %--------------------------------------------------------------------------    
SNMDA = [1   0   0   0   0;   %  SS
         0   1   0   0   0;   %  SP
         0   1   0   0   0;   %  SI
         0   0   0   0   0;   %  DP
         0   0   0   0   0;   %  DI
         0   0   0   0   0;   %  TP % 0 in ket study
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
GEa = [  0     0     0     0     0     2     0     2;
         2     2     0     0     0     0     0     0;
         0     2     0     0     0     0     0     0;
         0     2     0     0     0     0     0     0;
         0     0     0     2     0     0     0     0;
         0     0     0     2     0     0     0     0;
         0     0     0     0     0     0     0     2;
         2     0     0     0     0     2     0     0];

GEn =   [0     0     0     0     0     2     0     2;
         2     2     2     0     0     0     0     0;
         0     2     2     0     0     0     0     0;
         0     2     0     0     0     0     0     0;
         0     0     0     2     0     0     0     0;
         0     0     0     2     0     0     0     0;
         0     0     0     0     0     0     0     2;
         2     0     0     0     0     2     0     0];




% Inhibitory connections (np x np): GABA-A & GABA-B
%--------------------------------------------------------------------------
GIa =[8     0     10    0     0     0     0     0;
      0    18     10    0     0     0     0     0;
      0     0     10    0     0     0     0     0;
      0     0     0     8     6     0     0     0;
      0     0     0     0    14     0     0     0;
      0     0     0     0     6     8     0     0;
      0     0     0     0     0     0     8     0;
      0     0     0     0     0     0     8     8];

GIb = GIa;



% Channel rate constants [decay times]
%--------------------------------------------------------------------------
KE  = exp(-P.T(:,1))*1000/2.2;%3;            % excitatory rate constants (AMPA) % 2 to 5
KI  = exp(-P.T(:,2))*1000/5;%6;           % inhibitory rate constants (GABAa)
KN  = exp(-P.T(:,3))*1000/100;%40;          % excitatory rate constants (NMDA) 40-100
KB  = exp(-P.T(:,4))*1000/300;          % excitatory rate constants (NMDA)
KM  = (exp(-P.T(:,5))*1000/160) ;               % m-current opening + CV
KH  = (exp(-P.T(:,6))*1000/100) ;               % h-current opening + CV
Kht2a = (exp(-P.T(:,7))*1000/200) ;


% notes on time-constants:
%-----------------------------------------------------------------------
% cojuld even use number from this friston paper
%https://www.sciencedirect.com/science/article/pii/S0361923000004366?via%3Dihub
% ampa = 1.2 to 2.4 ms
% gabaa -   6ms
% nmda - 50 ms

% gaba-b maybe evern 300 or 500ms
% now using faster AMPA and GABA-A dynamics based on this book:
% https://neuronaldynamics.epfl.ch/online/Ch3.S1.html#:~:text=GABAA%20synapses%20have%20a,been%20deemed%203%20times%20larger.

% Trial-specific effects on time constants: AMPA & NMDA only for LTP task
if isfield(P,'T1')
    KE = KE + P.T1(1);
    KN = KN + P.T1(2);
end

% Voltages [reversal potentials] (mV)
%--------------------------------------------------------------------------
VL   = -70;                               % reversal  potential leak (K)
VE   =  60 ;                              % reversal  potential excite (Na)
VI   = -90 ;%* exp(P.pr(1));                % reversal  potential inhib (Cl)
VR   = -52 ;%* exp(P.pr(2));   %55          % threshold potential (firing)
VN   =  10 ;%* exp(P.pr(3));                % reversal Ca(NMDA)   
VB   = -100;%* exp(P.pr(4));               % reversal of GABA-B
VM   = -52;                            % reversal potential m-channels          
VH   = -30;                            % reversal potential h-channels 

% M- & H- channel conductances (np x np) {L6 & Thal Relay cells only}
%----------------------------------------------------------------------
% https://www.sciencedirect.com/science/article/pii/S0006349599769250
GIm = diag(4*[1 1 1 1 1 1 1 1].*exp(P.Mh(:)'));
GIh = diag(4*[0 0 0 0 0 1 0 1].*exp(P.Hh(:)'));

h     = 1 - spm_Ncdf_jdw(x(:,:,1),-100,300); % mean firing for h-currents
h     = 1 - 1./(1 + exp(-(2/3).*(x(:,:,1)-VH)));

Ght2a = diag([0 0 0 4 0 0 0 0]) * exp(P.HT2A);


% membrane capacitances {ss  sp  ii  dp  di  tp   rt  rl}
%--------------------------------------------------------------------------
CV   = exp(P.CV).*      [128*3 128 128/2 128 64  128  64  64*2]/1000;  

% leak conductance - fixed
%--------------------------------------------------------------------------
GL   = 1 ;       

% mean-field effects:
%==========================================================================

VR = VR + exp(P.S);

% neural-mass approximation to covariance of states: trial specific
%----------------------------------------------------------------------
R  = 2/3; %* exp(P.S); % P.S is the slope pf the sigmoid for each pop firing rate
FF = 1./(1 + exp(-R.*(x(:,:,1)-VR)));

RS = 30 ;
Fu = find( x(:,:,1) >= VR ); FF(Fu) = 1;
Fl = find( x(:,:,1) >= RS ); FF(Fl) = 0;
m  = FF;


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

% flow over every (ns x np) subpopulation
%==========================================================================
f     = x;

% Thalamo-cortical flow [eq. motion] over modes, populations, states...
%--------------------------------------------------------------------------
for i = 1:ns
               
        % input scaling: 
        %------------------------------------------------------------------
        dU = u(:)*C(i,1);
                                                
        % intrinsic coupling - parameterised
        %------------------------------------------------------------------
        E      = ( G(:,:,i).*GEa)*m(i,:)'; % AMPA currents
        ENMDA  = (Gn(:,:,i).*GEn)*m(i,:)'; % NMDA currents
        I      = ( G(:,:,i).*GIa)*m(i,:)'; % GABA-A currents
        IB     = (Gb(:,:,i).*GIb)*m(i,:)'; % GABA-B currents
                
            
        % intrinsic coupling - non-parameterised: intrinsic dynamics
        %--------------------------------------------------------------
        Im     = GIm*m(i,:)'; % M currents
        Ih     = GIh*h(i,:)'; % H currents

        Iht2a = Ght2a*m(i,:)';

        % extrinsic coupling (excitatory only) and background activity
        %------------------------------------------------------------------
        E     = (E     +  BE  + SA   *a (i,:)')*2;
        ENMDA = (ENMDA +  BE  + SNMDA*an(i,:)')*2;
                   
        % and exogenous input(U): 
        %------------------------------------------------------------------
        % flag for the oscillation injection desribed here: 
        % https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5310631/
        
        % if length of input vector > 1, represents more than one exogenous
        % input - one to thal relay and one to cortical pyramids
        if length(u) > 1
            E(8) = E(8) + dU(2);
            E(2) = E(2) + dU(1);
        else
            % otherwise just drive the thalamus
            input_cell        = [8 7];
            E(input_cell)     = E(input_cell) + dU;            
        end

        % direct current to thalamus
        % if isfield(P,'thi');
        %     E(8) = E(8) + exp(P.thi);
        %     ENMDA(8) = ENMDA(8) + exp(P.thi);
        % end

        % Scale NMDA conductance in L5 pyramidal cells based on 5HT-2A
        %x(i,4,4) = real( x(i,4,4) * (1 + 0.2 * x(i,4,8)) );  % 20% NMDA enhancement

        % --- 5HT-2A Effects on Other Conductances ---
        %if any(x(i,4,8) > 0)  % Only apply if serotonin is active
            % NMDA enhancement
            x(i,4,4) = x(i,4,4) * (1 + exp(P.s(1)) * 0.2 * x(i,4,8));

            % GABA-A suppression (disinhibition)
            x(i,3,3) = real( x(i,3,3) * (1 - exp(P.s(1)) * 0.2 * x(i,4,8)) );

            % GABA-B enhancement
            x(i,5,5) = x(i,5,5) * (1 + exp(P.s(1)) * 0.1 * x(i,4,8));

            % Kv7 (M-Channel) suppression (increased excitability)
            x(i,4,6) = x(i,4,6) * (1 - exp(P.s(1)) * 0.2 * x(i,4,8));

            % HCN (Ih) enhancement
            x(i,4,7) = x(i,4,7) * (1 + exp(P.s(1)) * 0.1 * x(i,4,8));
        %end
                              
        % Voltage equations
        %==================================================================
            
          % alternative magnesium block:
          %warning off;
          %mag_block = 1/(1 + 0.2*exp(-0.062*(exp(P.scale_NMDA))*squeeze(x(i,:,1))')) ;
          warning('off','all') ;
          mag_block = mldivide((1 + 0.2*exp(-0.062*(exp(P.scale_NMDA))*squeeze(x(i,:,1))'))',1)';
          %warning on;
          [~,warnID] = lastwarn;
          warning('off',warnID);

          f(i,:,1) =  (GL*(VL - x(i,:,1))+...
                       x(i,:,2).*((VE - x(i,:,1)))+...
                       x(i,:,3).*((VI - x(i,:,1)))+...
                       x(i,:,5).*((VB - x(i,:,1)))+...
                       x(i,:,6).*((VM - x(i,:,1)))+...
                       x(i,:,7).*((VH - x(i,:,1)))+...
                       x(i,:,4).*((VN - x(i,:,1))).*mag_block + ...
                       x(i,:,8).*((-10 - x(i,:,1))) )./CV;
          
                   
        % Conductance equations
        %==================================================================           
        
        f(i,:,2) = (E'     - x(i,:,2)).* (KE(i,:)');
        f(i,:,3) = (I'     - x(i,:,3)).* (KI(i,:)');
        f(i,:,5) = (IB'    - x(i,:,5)).* (KB(i,:)');
        f(i,:,4) = (ENMDA' - x(i,:,4)).* (KN(i,:)');        
        f(i,:,6) = (Im'    - x(i,:,6)).*(KM(i,:) );
        f(i,:,7) = (Ih'    - x(i,:,7)).*(KH(i,:) );
        f(i,:,8) = (Iht2a'  - x(i,:,8)).*(Kht2a');

end


% vectorise equations of motion
%==========================================================================
f = spm_vec((f));
pE = P;
 
[J,Q,D]=deal([]);

if (nargout < 2 || nargout == 50) && nargin < 5, return, end

% Only compute Jacobian (gradients) if requested
%==========================================================================
J = spm_cat(spm_diff(M.f,x,u,P,M,1));

if nargout < 3 && nargin < 5, return, end

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
d  = D.*full(exp(P.D(1:2)))/1000;
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

Tc              = zeros(np,np);
Tc([7 8],[1:6]) = CT  * exp(P.CT); % L6->thal
Tc([1:6],[7 8]) = TC  * exp(P.TC); % thal->ss

Tc = Tc / 1000;
Tc = kron(ones(nk,nk),kron(Tc,ones(ns,ns)));


%kd = exp(P.a(1)) * 8;
%ID = [4 1/4 1 8 1/2 4 2 20]/8;%2.4;
ID = [2 1 1 1 1 2 1 2];
ID = ID.*exp(P.ID)/1000; 
ID = repmat(ID,[1 nk]);

ID = repmat(ID(:)',[np*nk,1]);
ID = kron(ID,ones(ns,ns));

%ID = ID - ID(:);

% Mean intra-population delays, inc. axonal etc. Seem to help oscillation
%--------------------------------------------------------------------------
Dp = ~Ss;                            % states: different sources
Ds = ~Sp & Ss;                       % states: same source different pop.
%Ds = Ds.*(~(Ds & Tc));              % remove t-c and c-t from intrinsic

D = d(1)*Ds + Tc + (ID) ;
D =  Tc + (ID) ;

% Implement: dx(t)/dt = f(x(t - d)) = inv(1 - D.*dfdx)*f(x(t))
%                     = Q*f = Q*J*x(t)
%--------------------------------------------------------------------------
%Q  = spm_inv(speye(length(J)) - D.*J);
%Q  = spm_inv(D.*J);
