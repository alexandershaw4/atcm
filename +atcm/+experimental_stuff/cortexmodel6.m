function [f,J,Q] = cortexmodel6(x,u,P,M,m)
% State equations for an extended canonical thalamo-cortical neural-mass model.
%
% This model implements a conductance-based canonical cortical circuit,
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
% FORMAT [f,J,Q,D] = atcm.tcm_nomh(x,u,P,M)
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
%
%
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
Gn = full(P.Hn);
Gn = exp(Gn);

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
%       SP  DP  tp 
SA   = [1   0   0;   %  SS
        0   1   0;   %  SP
        0   1   0;   %  SI
        1   0   0;   %  DP
        0   0   0;   %  DI
        0   0   1 ]; %  TP % 
        
    
% % extrinsic NMDA-mediated connections (F B) - from superficial and deep pyramidal cells
% %--------------------------------------------------------------------------    
SNMDA =[1   0   0;   %  SS
        0   1   0;   %  SP
        0   1   0;   %  SI
        1   0   0;   %  DP
        0   0   0;   %  DI
        0   0   1 ]; %  TP 
     
% intrinsic connectivity switches
%--------------------------------------------------------------------------    
%   population: 1  - Spint stellates (L4)                : e
%               2  - Superficial pyramids (L2/3)         : e
%               3  - Inhibitory interneurons (L2/3)      : i
%               4  - Deep pyramidal cells (L5)           : e
%               5  - Deep interneurons (L5)              : i
%               6  - Thalamic projection neurons -L6     : e

GEa = zeros(6,6);
GIa = zeros(6,6);

% Excitatory (np x np): AMPA & NMDA
%--------------------------------------------------------------------------
% This is a simplified, predictive-coding friendly excitatory architecture
%           ss  sp  si  dp  di  tp    
GEa(1,:) = [0   0   0   0   0   2]/1;
GEa(2,:) = [4   0   0   0   0   0]/1; 
GEa(3,:) = [4   4   0   0   0   0]/1; 
GEa(4,:) = [0   4   0   0   0   0]/1;
GEa(5,:) = [0   0   0   4   0   0]/1;
GEa(6,:) = [0   0   0   2   0   0]/1; 


GEa = GEa .* ~eye(np);
%GEa = GEa * .8;
GEa = GEa + eye(np);      % KILLED
GEn = GEa;

GEn = GEn + (eye(np)/8);

% Inhibitory connections (np x np): GABA-A & GABA-B
%--------------------------------------------------------------------------
%           ss  sp  si  dp  di  tp  
GIa(1,:) = [8   0   8   0   0   0];
GIa(2,:) = [0   16  32  0   0   0];
GIa(3,:) = [0   0   32  0   0   0];
GIa(4,:) = [0   0   0   8   12  0];
GIa(5,:) = [0   0   0   0   16  0];
GIa(6,:) = [0   0   0   0   8   8];

GIb      = GIa;

if IncludeMH
    
    % M- & H- channel conductances (np x np) {L6 & Thal Relay cells only}
    %----------------------------------------------------------------------
    VM   = -70;                            % reversal potential m-channels          
    VH   = -30;                            % reversal potential h-channels 

    GIm  = diag(6)/4;
    Mh   = diag(exp(P.Mh));
    GIh  = diag(6)/4;
    GIh  = diag(exp(P.Hh));

    KM   = (exp(-P.m)*1000/160) ;               % m-current opening + CV
    KH   = (exp(-P.h)*1000/100) ;               % h-current opening + CV
    h    = 1 - spm_Ncdf_jdw(x(:,:,1),-100,300); % mean firing for h-currents
end

% Channel rate constants [decay times]
%--------------------------------------------------------------------------
KE  = exp(-P.T(:,1))*1000/4;            % excitatory rate constants (AMPA)
KI  = exp(-P.T(:,2))*1000/16;           % inhibitory rate constants (GABAa)
KN  = exp(-P.T(:,3))*1000/100;          % excitatory rate constants (NMDA)
KB  = exp(-P.T(:,4))*1000/200;          % excitatory rate constants (NMDA)

%KE  = exp(-P.T(:,1))*1000/5;            % excitatory rate constants (AMPA)
KN  = exp(-P.T(:,3))*1000/150;          % excitatory rate constants (NMDA)
%KI  = exp(-P.T(:,2))*1000/4;           % inhibitory rate constants (GABAa)

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
VR   = -40;                               % threshold potential (firing)
VN   =  10;                               % reversal Ca(NMDA)   
VB   = -100;                              % reversal of GABA-B

% membrane capacitances {ss  sp  ii  dp  di  tp}
%--------------------------------------------------------------------------
CV   = exp(P.CV).*      [128 128 64  128 64  128]/1000;  

% leak conductance - fixed
%--------------------------------------------------------------------------
GL   = 1;          

% mean-field effects:
%==========================================================================

% neural-mass approximation to covariance of states: trial specific
%----------------------------------------------------------------------
Vx   = exp(P.S)*32;
if nargin < 5
    % compute only if not passed by integrator
    m    =     spm_Ncdf_jdw(x(:,:,1),VR,Vx);
end

% extrinsic effects
%--------------------------------------------------------------------------
a       = zeros(ns,3);
an      = zeros(ns,3); 

a(:,1)  = A{1}*m(:,2);                      % forward afference  AMPA - SP->SS
a(:,2)  = A{2}*m(:,4);                      % backward afference AMPA - DP->SP
a(:,3)  = A{3}*m(:,6);                      % FWD thalamic projection pyramids
an(:,1) = AN{1}*m(:,2);                     % forward afference  NMDA
an(:,2) = AN{2}*m(:,4);                     % backward afference NMDA
an(:,3) = AN{3}*m(:,6);                     % thalamic projection pyramids

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
        if any(full(U(:))) ; %&& size(U,1) >= i
            %dU = u(1)*( C(i,:).*[1 1/64 1/128 1/128] );                  % CHANGE BACK
            dU = u*C(i,1);%1
        else
            dU = 0;
        end
        
        % intrinsic coupling - parameterised
        %------------------------------------------------------------------
        E      = ( G(:,:,i).*GEa)*m(i,:)'; % AMPA currents
        ENMDA  = ( Gn(:,:,i).*GEn)*m(i,:)'; % NMDA currents
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
        E     = (E     +  BE  + SA*a(i,:)')*2;
        ENMDA = (ENMDA +  BE  + SNMDA*an(i,:)')*2;
      
        % and exogenous input(U): 
        %------------------------------------------------------------------
        input_cell        = 1;
                
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
        %pop_rates = [1 1 2 1 2 1 1 1];
        %pop_rates = pop_rates.*exp(P.pr);
        
        f(i,:,2) = (E'     - x(i,:,2)).* KE(i,:);%*pop_rates);
        f(i,:,3) = (I'     - x(i,:,3)).* KI(i,:);%*pop_rates);
        f(i,:,5) = (IB'    - x(i,:,5)).* KB(i,:);%*pop_rates);
        f(i,:,4) = (ENMDA' - x(i,:,4)).* KN(i,:);%*pop_rates);
        
        if IncludeMH
            f(i,:,6) = (Im'    - x(i,:,6)).*(KM(i,:) );
            f(i,:,7) = (Ih'    - x(i,:,7)).*(KH(i,:) );
        end

        % Restrict the states flow rate by the conduction rate of the population
        % c.f. synaptic delays + conduction delays
        %------------------------------------------------------------------
        DV       = 1./[1 1 1 1 1 1]; 
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
D  = [.6 16];
d  = -D.*full(exp(P.D(1:2)))/1000;
Sp = kron(ones(nk,nk),kron( eye(np,np),eye(ns,ns)));  % states: same pop.
Ss = kron(ones(nk,nk),kron(ones(np,np),eye(ns,ns)));  % states: same source


if isfield(P,'ID')
    % ignore..... doesn't trigger  unless you have an entry in P 'ID'
    %-----------------------------------------------------------------
    % intrisc delays
    ID = [0 0 0 1 0 1 1 1];
    ID = [1 .2 .1 1 .2 1 .4 1];
    ID = [2 1  .1 2 .2 2]; % this 
    ID = -ID.*exp(P.ID)/1000;
    ID = kron(ones(nk,nk),kron(diag(ID),eye(ns,ns)));
    Tc = ID;
end


% Mean intra-population delays, inc. axonal etc. Seem to help oscillation
%--------------------------------------------------------------------------
Dp = ~Ss;                            % states: different sources
Ds = ~Sp & Ss;                       % states: same source different pop.
%Ds = Ds.*(~(Ds & Tc));              % remove t-c and c-t from intrinsic
D  = d(2)*Dp + d(1)*Ds + Tc  ;       %+ Dself;% Complete delay matrix

% Implement: dx(t)/dt = f(x(t - d)) = inv(1 - D.*dfdx)*f(x(t))
%                     = Q*f = Q*J*x(t)
%--------------------------------------------------------------------------
Q  = spm_inv(speye(length(J)) - D.*J);

