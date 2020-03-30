function [f,J,Q,D] = tcm_traub(x,u,P,M)
% State equations for an extended canonical thalamo-cortical neural-mass model.
%
% This model implements a conductance-based canonical thalamo-cortical circuit,
% with cytoarchitecture inspired by Gilbert & Wiesel (1983), Douglas & 
% Martin (2004) and Traub (2004) models.
%
% The equations of motion are Moris Lecar-esque equations, similar to Moran
% (2011), but with conductances for AMPA, NMDA, GABA-A, GABA-B, m- and h-
% channels. These 'channels' feature their own reversal poentials and rate
% constants:
%
% K  = -70           (Leak)
% Na =  60  / 4 ms   (AMPA)
% Cl = -90  / 16 ms  (GABA-A)
% Ca =  10  / 100 ms (NMDA)   + voltage mag switch
% B  = -100 / 200 ms (GABA-B)
% m  = -70  / 160 ms (m-chan)
% h  = -30  / 100 ms (h-chan)
% firing = -40
%
% FORMAT [f,J,Q,D] = spm_fx_cmm_NMDA_ExtendedV1mh(x,u,P,M)
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
%        state: 1 V  - voltage
%               2 gE - conductance: AMPA   (excitatory)
%               3 gI - conductance: GABA-A (inhibitory)
%               4 gN - conductance: NMDA   (excitatory)
%               5 gB - conductance: GABA-B (inhibitory)
%               6 gm - conductance: m-channels (noninactivating)
%               7 gh - conductance: h-channels
%
%      outputs: f = model states as a vector - hint: spm_unvec(f,M.x) 
%               J = system Jacobian - dfdx
%               Q = delay operator  - Q = inv(1 - D.*dfdx)*f(x(t))
%               D = states delay matrix
%
%
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
    if i > 2
        A{i}  = A{i}  *0;
        AN{i} = AN{i} *0;
    end
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


%--------------------------------------------------------------------------
G    = full(P.H);
G    = exp(G);


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
SA   = [1   0   0   0   0;   %  SS
        0   1   0   0   0;   %  SP
        0   1   0   0   0;   %  SI
        1   0   0   0   0;   %  DP
        0   0   0   0   0;   %  DI
        0   0   0   0   0;   %  TP
        0   0   0   0   0;   %  rt
        0   0   0   0   0]/8;%  rc
    
% % extrinsic NMDA-mediated connections (F B) - from superficial and deep pyramidal cells
% %--------------------------------------------------------------------------    
SNMDA = [1   0   0   0   0;   %  SS
         0   1   0   0   0;   %  SP
         0   1   0   0   0;   %  SI
         1   0   0   0   0;   %  DP
         0   0   0   0   0;   %  DI
         0   0   0   0   0;   %  TP
         0   0   0   0   0;   %  rt
         0   0   0   0   0]/8;%  rc

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


% Excitatory (np x np): AMPA & NMDA
%--------------------------------------------------------------------------
%           ss  sp  si  dp  di  tp  rt  rl   
% GEa(1,:) = [0   .1  0   .5  0   0   0   1]/1;% ss
% GEa(2,:) = [1   0   0   .5  0   .5  0   .8]/1;% sp
% GEa(3,:) = [1   3   0   1   0   0   0   1]/1;% si
% GEa(4,:) = [1   .1  0   0   0   0   0   1.5]/1;% dp
% GEa(5,:) = [1   1   0   3   0   0   0   1.5]/1;% di
% GEa(6,:) = [0   .5  0   2   0   0   0   1]/1;% tp
% GEa(7,:) = [0   0   0   0   0   .5  0   0]/1;% rt
% GEa(8,:) = [0   0   0   0   0   .75 0   0]/1;% rl

GEa(1,:) = [0   4  0   4   0   0   0   4]/1;% ss
GEa(2,:) = [4   0  0   4   0   4   0   4]/1;% sp
GEa(3,:) = [4   4  0   4   0   0   0   4]/1;% si
GEa(4,:) = [4   4  0   0   0   0   0   4]/1;% dp
GEa(5,:) = [4   4  0   4   0   0   0   4]/1;% di
GEa(6,:) = [0   4  0   4   0   0   0   4]/1;% tp
GEa(7,:) = [0   0  0   0   0   4   0   0]/1;% rt
GEa(8,:) = [0   0  0   0   0   4   0   0]/1;% rl


GEa = GEa .* ~eye(np);
GEa = GEa + eye(np);

% Append self-excitation parameters {excitatory cells only}
%--------------------------------------------------------------------------
for i = 1:ns
    sG(:,i) = exp(P.G(i,:)) .* ([4 4 0 4 0 4 0 4]);
end

% NMDA conductances: same as AMPA
%--------------------------------------------------------------------------
GEn = GEa;

   
% Inhibitory connections (np x np): GABA-A & GABA-B
%--------------------------------------------------------------------------
%           ss  sp  si  dp  di  tp  rt  rl
% GIa(1,:) = [8   0   .1  0   1.5 0   0   0 ];% ss
% GIa(2,:) = [0   8   1.2 0   0   0   0   0 ];% sp
% GIa(3,:) = [0   0   .2  0   0   0   0   0 ];% si
% GIa(4,:) = [0   0   0   8   1   0   0   0 ];% dp
% GIa(5,:) = [0   0   0   0   .2  0   0   0 ];% di
% GIa(6,:) = [0   0   0   0   .7  8   0   0 ];% tp
% GIa(7,:) = [0   0   0   0   0   0   8   0 ];% rt
% GIa(8,:) = [0   0   0   0   0   0   1.5 8 ];% rl

GIa(1,:) = [8   0   2   0   8   0   0   0 ];% ss
GIa(2,:) = [0   8   2   0   0   0   0   0 ];% sp
GIa(3,:) = [0   0  32   0   0   0   0   0 ];% si
GIa(4,:) = [0   0   0   8   8   0   0   0 ];% dp
GIa(5,:) = [0   0   0   0  32   0   0   0 ];% di
GIa(6,:) = [0   0   0   0   8   8   0   0 ];% tp
GIa(7,:) = [0   0   0   0   0   0  32   0 ];% rt
GIa(8,:) = [0   0   0   0   0   0   8   8 ];% rl


% GABA-B conductances: same as GABA-A
%--------------------------------------------------------------------------
GIb     = GIa;

% M- channel connections (np x np) {excitatory cells only}
%--------------------------------------------------------------------------
GIm  = eye(8)*4/10;
Mh   = diag(exp(P.Mh));
 
% H-channel conductances {L6 & Thal Relay only}
%--------------------------------------------------------------------------
GIh      = full(sparse([6 8],[6 8],4/10,8,8));
Hh       = exp(P.Hh);
GIh(6,6) = GIh(6,6)*Hh(1);
GIh(8,8) = GIh(8,8)*Hh(2);

% Channel rate constants [opening times] & conduction velocities (ns x np) 
%--------------------------------------------------------------------------
RR = exp(-P.T) .* (1000./[4 16 100 200]);           % receptor opening rate
CR = exp(-P.TV).* (100./[10 2 2 10 5 10 10 800]);
CR = CR*0;

for i = 1:8
    for j = 1:4
        Kmat(i,j) = RR(j) + CR(i);
    end
end

for in = 1:ns
    KE(in,:)    = Kmat(:,1)';        % AMPA   opening + CV
    KI(in,:)    = Kmat(:,2)';        % GABA-A opening + CV
    KN(in,:)    = Kmat(:,3)';        % NMDA   opening + CV
    KB(in,:)    = Kmat(:,4)';        % GABA-B opening + CV
end


% % AMPA rates
% %     ss sp si dp di tp rt rl
% KEt = [0  2  0  2  0  2  0  2; % ss
%       2  0  0  2  0  2  0  2; % sp
%       .8 .8 0  .8 0  .8 0  1; % si
%       2  2  0  0  0  2  0  2; % dp
%       .8 .8 0  .8 0  .8 0  1; % di
%       2  2  0  2  0  0  0  2; % tp
%       0  0  0  0  0  0  0  2; % rt
%       0  0  0  0  0  2  3.3 0];% rl
%   
% % NMDA rates
% %     ss sp si dp di tp rt rl
% KNt = [0  13 0  13 0  13 0  13;
%       13 0  0  13 0  13 0  13;
%       0  10 0  10 0  10 0  10;
%       13 13 0  0  0  13 0  13;
%       0  10 0  10 0  10 0  10;
%       13 13 0  13 0  0  0  13;
%       0  0  0  0  0  0  0  15;
%       0  0  0  0  0  0  0  0]*10;
% 
% % GABA-A rates
% %     ss sp si dp di tp rt rl
% KIt = [0  0  3  0  3  0  0  0;
%       0  0  3  0  3  0  0  0;
%       0  0  0  0  0  0  0  0;
%       0  0  3  0  3  0  0  0;
%       0  0  0  0  0  0  0  0;
%       0  0  3  0  3  0  0  0;
%       0  0  0  0  0  0  9  0;
%       0  0  0  0  0  0  3.3 0];
% 
% % invert & multiply by parameter
% for in = 1:ns
%     KE(in,:,:) = exp( -P.T(1) ) * (1000./KEt);
%     KI(in,:,:) = exp( -P.T(2) ) * (1000./KIt);
%     KN(in,:,:) = exp( -P.T(3) ) * (1000./KNt);
%     KB(in,:,:) = exp( -P.T(4) ) * (1000./KIt);
% end
% 
% KE(isinf(KE)) = 0;
% KI(isinf(KI)) = 0;
% KN(isinf(KN)) = 0;
% KB(isinf(KB)) = 0;


KM    = (exp(-P.m)*1000/160) ;   % m-current opening + CV
KH    = (exp(-P.h)*1000/100) ;   % h-current opening + CV

% Voltages [reversal potentials] (mV)
%--------------------------------------------------------------------------
VL   = -70;                               % reversal  potential leak (K)
VE   =  60;                               % reversal  potential excite (Na)
VI   = -90;                               % reversal  potential inhib (Cl)
VR   = -40;                               % threshold potential (firing)
VN   =  10;                               % reversal Ca(NMDA)   
VB   = -100;                              % reversal of GABA-B
VM   = -70;                               % reversal potential m-channels          
VH   = -30;                               % reversal potential h-channels 

% membrane capacitances {ss  sp  ii  dp  di  tp   rt  rl}
%--------------------------------------------------------------------------
CV   = exp(P.CV).*[128 128 128 128 128 128 128 128]/1000;  
CV   = exp(P.CV).*[128 128 256 128 256 128 256 128]/1000;  
CV   = exp(P.CV).*[16 16 32 16 32 16 32 16]/1000;

% leak conductance - fixed
%--------------------------------------------------------------------------
GL   = 1;          

% mean-field effects:
%==========================================================================

% neural-mass approximation to covariance of states: trial specific
%----------------------------------------------------------------------
Vx   = exp(P.S)*32;
h    = 1 - spm_Ncdf_jdw(x(:,:,1),-100,300); % mean firing for h-currents
m    =     spm_Ncdf_jdw(x(:,:,1),VR,Vx);    % mean firing rate  

% extrinsic effects
%--------------------------------------------------------------------------
a(:,1)  = A{1}*m(:,2);                      % forward afference  AMPA
a(:,2)  = A{2}*m(:,4);                      % backward afference AMPA
a(:,3)  = A{3}*m(:,6);                      % thalamic projection pyramids
a(:,4)  = A{4}*m(:,7);                      % reticular AMPA
a(:,5)  = A{5}*m(:,8);                      % relay AMPA
an(:,1) = AN{1}*m(:,2);                     % forward afference  NMDA
an(:,2) = AN{2}*m(:,4);                     % backward afference NMDA
an(:,3) = AN{3}*m(:,6);                     % thalamic projection pyramids
an(:,4) = AN{4}*m(:,7);                     % reticular NMDA
an(:,5) = AN{5}*m(:,8);                     % relay NMDA

% Averge background activity and exogenous input
%==========================================================================
BE     = exp(P.E)*0.8;

% input
%--------------------------------------------------------------------------
if isfield(M,'u')
      U =   u(:); % endogenous input
else; U = C*u(:); % exogenous input
end

% flow over every (ns x np) subpopulation
%==========================================================================
f     = x;

% Thalamo-cortical flow [motion] over modes, populations, states...
%--------------------------------------------------------------------------
for i = 1:ns
   
        % intrinsic coupling
        %------------------------------------------------------------------
        EG = ( ~eye(np) .* G(:,:,i) ) + diag(sG(:,i));
        
        E      = (EG       .*GEa)*m(i,:)'; % AMPA currents
        ENMDA  = (EG       .*GEn)*m(i,:)'; % NMDA currents
        I      = ( G(:,:,i).*GIa)*m(i,:)'; % GABA-A currents
        IB     = ( G(:,:,i).*GIb)*m(i,:)'; % GABA-B currents
        Im     = (Mh(:,:,i).*GIm)*m(i,:)'; % m currents
        Ih     =             GIh *h(i,:)'; % h currents
        
        % extrinsic coupling (excitatory only) and background activity
        %------------------------------------------------------------------
        E     = (E     +  BE  + SA   *a (i,:)')*2;
        ENMDA = (ENMDA +  BE  + SNMDA*an(i,:)')*2;
      
        % and exogenous input(U): 
        %------------------------------------------------------------------
        input_cell = 8; 
        
        try   E(input_cell)          = E(input_cell)           + U(i,1);
              ENMDA(input_cell)      = ENMDA(input_cell)       + U(i,1);
        catch E(input_cell)          = E(input_cell)           + U(1,1);
              ENMDA(input_cell)      = ENMDA(input_cell)       + U(1,1);
        end
        
        % Voltage equation
        %==================================================================
          f(i,:,1) =         (GL*(VL - x(i,:,1))+...
                       x(i,:,2).*(VE - x(i,:,1))+...
                       x(i,:,3).*(VI - x(i,:,1))+...
                       x(i,:,6).*(VM - x(i,:,1))+...
                       x(i,:,7).*(VH - x(i,:,1))+...
                       x(i,:,5).*(VB - x(i,:,1))+...
                       x(i,:,4).*(VN - x(i,:,1)).*mg_switch(x(i,:,1)))./CV;
                   
        % Conductance equation
        %==================================================================        
        f(i,:,2) = (E'     - x(i,:,2)).*KE(i,:);
        f(i,:,3) = (I'     - x(i,:,3)).*KI(i,:);
        f(i,:,5) = (IB'    - x(i,:,5)).*KB(i,:);
        f(i,:,4) = (ENMDA' - x(i,:,4)).*KN(i,:);
        f(i,:,6) = (Im'    - x(i,:,6)).*(KM(i,:) );
        f(i,:,7) = (Ih'    - x(i,:,7)).*(KH(i,:) );
        
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


% Mean intra-population delays, inc. axonal etc. Seem to help oscillation
%--------------------------------------------------------------------------
Dp = ~Ss;                            % states: different sources
Ds = ~Sp & Ss;                       % states: same source different pop.
%Ds = Ds.*(~(Ds & Tc));               % remove t-c and c-t from intrinsic
D  = d(2)*Dp + d(1)*Ds;       %+ Dself;% Complete delay matrix

% Implement: dx(t)/dt = f(x(t - d)) = inv(1 - D.*dfdx)*f(x(t))
%                     = Q*f = Q*J*x(t)
%--------------------------------------------------------------------------
Q  = spm_inv(speye(length(J)) - D.*J);


