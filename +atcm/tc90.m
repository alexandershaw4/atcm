function [f,J,Q,D] = tc90(x,u,P,M,m)
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
%
%
%--------------------------------------------------------------------------

 
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


% % % extrinsic connections (F B) - from superficial and deep pyramidal cells
% % %--------------------------------------------------------------------------
% %       SP  DP  tp  rt  rc
% SA   = [1   0   1   0   0;   %  SS    % added TP->SP
%         0   1   0   0   0;   %  SP
%         0   1   0   0   0;   %  SI
%         1   0   0   0   0;   %  DP
%         0   0   0   0   0;   %  DI
%         0   0   1   0   0;   %  TP % 0 in ket study
%         0   0   0   0   1;   %  rt % 0 in ket study
%         0   0   0   1   0]/8;%  rc % 0 in ket study
%     
%     %SA(:,[3 4 5]) = 0; % For ket study
%     
%     
% % % extrinsic NMDA-mediated connections (F B) - from superficial and deep pyramidal cells
% % %--------------------------------------------------------------------------    
% SNMDA = [1   0   1   0   0;   %  SS
%          0   1   0   0   0;   %  SP
%          0   1   0   0   0;   %  SI
%          1   0   0   0   0;   %  DP
%          0   0   0   0   0;   %  DI
%          0   0   1   0   0;   %  TP % 0 in ket study
%          0   0   0   0   1;   %  rt % 0 in ket study
%          0   0   0   1   0]/8;%  rc % 0 in ket study
% 
%      %SNMDA(:,[3 4 5]) = 0; % For ket study
     
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

GEa = zeros(90,90);
GIa = zeros(90,90);

[M,C,L,B] = AALMODELS;

GEa = M{1}.F | L; % ampa fwd
GEn = M{1}.B; % nmda bkw

C = C(:).*exp(P.C(:));

GIa = GEa;
GIb = GIa;


% Channel rate constants [decay times]
%--------------------------------------------------------------------------
KE  = exp(-P.T(:,1))*1000/4;            % excitatory rate constants (AMPA)
KI  = exp(-P.T(:,2))*1000/16;           % inhibitory rate constants (GABAa)
KN  = exp(-P.T(:,3))*1000/100;          % excitatory rate constants (NMDA)
KB  = exp(-P.T(:,4))*1000/200;          % excitatory rate constants (NMDA)

% KE  = exp(-P.T(:,:,1))*1000/4;            % excitatory rate constants (AMPA)
% KI  = exp(-P.T(:,:,2))*1000/16;           % inhibitory rate constants (GABAa)
% KN  = exp(-P.T(:,:,3))*1000/100;          % excitatory rate constants (NMDA)
% KB  = exp(-P.T(:,:,4))*1000/200;          % excitatory rate constants (NMDA)

% now using faster AMPA and GABA-A dynamics based on this book:
% https://neuronaldynamics.epfl.ch/online/Ch3.S1.html#:~:text=GABAA%20synapses%20have%20a,been%20deemed%203%20times%20larger.

%KE  = exp(-P.T(:,1))*1000/3;            % excitatory rate constants (AMPA)
%KN  = exp(-P.T(:,3))*1000/150;          % excitatory rate constants (NMDA)
KI  = exp(-P.T(:,2))*1000/6;           % inhibitory rate constants (GABAa)


% Voltages [reversal potentials] (mV)
%--------------------------------------------------------------------------
VL   = -70;                               % reversal  potential leak (K)
VE   =  60;                               % reversal  potential excite (Na)
VI   = -90;                               % reversal  potential inhib (Cl)
VR   = -40;%-40;                               % threshold potential (firing)
VN   =  10;                               % reversal Ca(NMDA)   
VB   = -100;                              % reversal of GABA-B

% membrane capacitances {ss  sp  ii  dp  di  tp   rt  rl}
%--------------------------------------------------------------------------
CV   = exp(P.CV).*      128*ones(1,90)/1000;  


% leak conductance - fixed
%--------------------------------------------------------------------------
GL   = 1;          

% mean-field effects:
%==========================================================================

% neural-mass approximation to covariance of states: trial specific
%----------------------------------------------------------------------
Vx   = exp(P.S)*32; % 32
if nargin < 5
    % compute only if not passed by integrator
    m    =     spm_Ncdf_jdw(x(:,:,1),VR,Vx);
end


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
for i = 1;
            
        % input scaling: 
        %------------------------------------------------------------------
        if any(full(U(:))) ;
            dU = u(1)*C(i,1) ;
        else
            dU = 0;
        end
        
        %Gsc = ~eye(8);
        %Gsc = Gsc + (diag(exp(P.Gsc)));
        
        % intrinsic coupling - parameterised
        %------------------------------------------------------------------
        E      = ( G(:,:,i).*GEa)*m(i,:)'; % AMPA currents
        ENMDA  = (Gn(:,:,i).*GEn)*m(i,:)'; % NMDA currents
        I      = ( G(:,:,i).*GIa)*m(i,:)'; % GABA-A currents
        IB     = ( G(:,:,i).*GIb)*m(i,:)'; % GABA-B currents
                
        % extrinsic coupling (excitatory only) and background activity
        %------------------------------------------------------------------
        E     = (E     +  BE  )*2;
        ENMDA = (ENMDA +  BE  )*2;
      
        % and exogenous input(U): 
        %------------------------------------------------------------------
        input_cell        = 1:90;
        
        if isfield(M,'inputcell');
            input_cell = M.inputcell;
        end
        
        E(input_cell)     = E(input_cell)         +dU';
        ENMDA(input_cell) = ENMDA(input_cell)     +dU';
                
        % Voltage equation
        %==================================================================            
          f(i,:,1) =         (GL*(VL - x(i,:,1))+...
                       1.0*x(i,:,2).*(VE - x(i,:,1))+...
                       1.0*x(i,:,3).*(VI - x(i,:,1))+...
                       1.0*x(i,:,5).*(VB - x(i,:,1))+...
                       1.0*x(i,:,4).*(VN - x(i,:,1)).*mg_switch(x(i,:,1)))./CV;
            
                   
        % Conductance equations
        %==================================================================   
        %pop_rates = [1 1 1 1 1 1 1 1];
        %pop_rates = pop_rates.*exp(P.pr);
        
        %gabaa_rate = pop_rates .* exp(P.gaba);
        
        pop_rates = ones(1,90);
        gabaa_rate = pop_rates;
        
        f(i,:,2) = (E'     - x(i,:,2)).* (KE(i,:)*pop_rates);
        f(i,:,3) = (I'     - x(i,:,3)).* (KI(i,:)*gabaa_rate);
        f(i,:,5) = (IB'    - x(i,:,5)).* (KB(i,:)*pop_rates);
        f(i,:,4) = (ENMDA' - x(i,:,4)).* (KN(i,:)*pop_rates);
                
        
        % c.f. synaptic delays + conduction delays
        %------------------------------------------------------------------
        %DV       = 1./[1 1 1 2.2 1 2 8 8]; 
        %DV       = 1./[2 1 1 2.2 1 2 1 2]; 
        %DV       = 1./[1 1 2 1   2 1 1 1]; 
        
        %DV       = 1./[1 1 .2 2 .4 2 .8 1];
%         DV = 1./[1 1 1 1 1 1 1 1];
%         if isfield(P,'TV')
%             DV       = DV.*exp(P.TV);
%             f(i,:,2) = f(i,:,2) .* DV;  % AMPA
%             f(i,:,3) = f(i,:,3) .* DV;  % GABA-A
%             f(i,:,4) = f(i,:,4) .* DV;  % NMDA
%             f(i,:,5) = f(i,:,5) .* DV;  % GABA-B
% 
%         end
        
                
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

% if isfield(P,'ID')
%     % ignore..... doesn't trigger  unless you have an entry in P 'ID'
%     %-----------------------------------------------------------------
%     % intrisc delays
%     ID = [0 0 0 1 0 1 1 1];
%     ID = [1 .2 .1 1 .2 1 .4 1];
%     ID = [2 1  .1 2 .2 2 .4 2]; % this 
%        
%     ID = [1 1 .5 1 .5 1 .5 1]*(.6); % this 
%     
%     %ID = double(~~GEa | ~~GIa);
%     %ID = (repmat(ID,[8 1]).*~eye(8)+diag(ID)).* double(~~GEa | ~~GIa);
%     
%     %ID = (repmat(ID,[8 1])).* double(~~GEa | ~~GIa);
%     
%     %ID = diag(ID) + 1e-2*double(~~GEa | ~~GIa);
%     
%     ID = -ID.*exp(P.ID)/1000;
%     %ID = kron(ones(nk,nk),kron(diag(ID),eye(ns,ns)));
%     %IDm = ID+ID';
%     %IDm = IDm.*~eye(8);
%     %IDm = IDm + diag(ID);
%     IDm=ID;
%     ID = kron(ones(nk,nk),kron(diag(ID),eye(ns,ns)));
%     Tc = Tc + ID;
% end


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