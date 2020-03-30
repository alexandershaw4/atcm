function [f,J,Q,D,DV] = tcm_0(x,u,P,M,m)
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
%        state: 1 V  - voltage
%               2 gE - conductance: AMPA   (excitatory)
%               3 gI - conductance: GABA-A (inhibitory)
%               4 gN - conductance: NMDA   (excitatory)
%               5 gB - conductance: GABA-B (inhibitory)
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

% Flag: include M- & H- channels on L6 TP & Thalamic Relay cells, or not
%--------------------------------------------------------------------------
IncludeMH = 0;
 
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


% Channel rate constants [opening times] & conduction velocities 
%--------------------------------------------------------------------------
KE  = exp(-P.T(:,1))*1000/4;            % excitatory rate constants (AMPA)
KI  = exp(-P.T(:,2))*1000/16;           % inhibitory rate constants (GABAa)
KN  = exp(-P.T(:,3))*1000/100;          % excitatory rate constants (NMDA)
KB  = exp(-P.T(:,4))*1000/200;          % excitatory rate constants (NMDA)

% Voltages [reversal potentials] (mV)
%--------------------------------------------------------------------------
VL   = -70;                               % reversal  potential leak (K)
VE   =  60;                               % reversal  potential excite (Na)
VI   = -90;                               % reversal  potential inhib (Cl)
VR   = -40;                               % threshold potential (firing)
VN   =  10;                               % reversal Ca(NMDA)   
VB   = -100;                              % reversal of GABA-B

% membrane capacitances {ss  sp  ii  dp  di  tp   rt  rl}
%--------------------------------------------------------------------------
CV   = exp(P.CV).*      [128 32  32  128 64  128  64 32]/1000;  

% leak conductance - fixed
%--------------------------------------------------------------------------
GL   = 1;          

% neural-mass approximation to covariance of states: trial specific
%----------------------------------------------------------------------
Vx   = exp(P.S)*32;
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
      U =   u(:); % endogenous input
else; U = C*u(:); % exogenous input
end


% Delays
DV = 1./[1 1 1 1 1 1 8 8];
DV = 1./[1 1 1 2.2 1 2 8 8];
DV = DV.*exp(P.TV);


% motion
%==========================================================================
f = x;
i = 1;

if any( U(i,1) )
    U = C.*[1 1/64 1/128 1/128]*u(:);
else
    U = [0 0 0 0];
end

% L4 stellates
%==========================================================================
% L4: voltage equations
f(i,1,1) = (GL*(VL - x(i,1,1))+...
     x(i,1,2).*(VE - x(i,1,1))+...
     x(i,1,3).*(VI - x(i,1,1))+...
     x(i,1,5).*(VB - x(i,1,1))+...
     x(i,1,4).*(VN - x(i,1,1)).*mg_switch(x(i,1,1)))./CV(1);

% L4: spiny stellate: current equations
E = exp(P.H(1,:)).*[0   0   0   0   0   2   0   4] .* m(i,:);
I = exp(P.H(1,:)).*[8   0   2   0   0   0   0   0] .* m(i,:);
E = E + U(2) + BE;

E = sum(E);
I = sum(I);

% L4: spiny stellate: current equations
f(i,1,2) = ((E' - x(i,1,2)).*KE(i,:)) .* DV(1);
f(i,1,3) = ((I' - x(i,1,3)).*KI(i,:)) .* DV(1);
f(i,1,4) = ((E' - x(i,1,4)).*KN(i,:)) .* DV(1);
f(i,1,5) = ((I' - x(i,1,5)).*KB(i,:)) .* DV(1);


 
% L2/3 pyramids
%==========================================================================
%  L2/3 pyramids voltage equations
f(i,2,1) = (GL*(VL - x(i,2,1))+...
     x(i,2,2).*(VE - x(i,2,1))+...
     x(i,2,3).*(VI - x(i,2,1))+...
     x(i,2,5).*(VB - x(i,2,1))+...
     x(i,2,4).*(VN - x(i,2,1)).*mg_switch(x(i,2,1)))./CV(2);

%  L2/3 pyramids: current equations
E = exp(P.H(2,:)).*[4   0   0   0   0   0   0   0] .* m(i,:);
I = exp(P.H(2,:)).*[0   16  16  0   0   0   0   0] .* m(i,:);

E = sum(E) + U(3) + BE;
I = sum(I);

%  L2/3 pyramids: current equations
f(i,2,2) = ((E' - x(i,2,2)).*KE(i,:)) .* DV(2);
f(i,2,3) = ((I' - x(i,2,3)).*KI(i,:)) .* DV(2);
f(i,2,4) = ((E' - x(i,2,4)).*KN(i,:)) .* DV(2);
f(i,2,5) = ((I' - x(i,2,5)).*KB(i,:)) .* DV(2);


% L2/3 interneurons
%==========================================================================
%  L2/3 interneurons voltage equations
f(i,3,1) = (GL*(VL - x(i,3,1))+...
     x(i,3,2).*(VE - x(i,3,1))+...
     x(i,3,3).*(VI - x(i,3,1))+...
     x(i,3,5).*(VB - x(i,3,1))+...
     x(i,3,4).*(VN - x(i,3,1)).*mg_switch(x(i,3,1)))./CV(3);

%  L2/3 interneurons: current equations
E = exp(P.H(3,:)).*[4   4   0   0   0   0   0   0] .* m(i,:);
I = exp(P.H(3,:)).*[0   0   32  0   0   0   0   0] .* m(i,:);

E = sum(E) + BE;
I = sum(I);

%  L2/3 interneurons: current equations
f(i,3,2) = ((E' - x(i,3,2)).*KE(i,:)) .* DV(3);
f(i,3,3) = ((I' - x(i,3,3)).*KI(i,:)) .* DV(3);
f(i,3,4) = ((E' - x(i,3,4)).*KN(i,:)) .* DV(3);
f(i,3,5) = ((I' - x(i,3,5)).*KB(i,:)) .* DV(3);


% L5 pyramids
%==========================================================================
%  L5 pyramids voltage equations
f(i,4,1) = (GL*(VL - x(i,4,1))+...
     x(i,4,2).*(VE - x(i,4,1))+...
     x(i,4,3).*(VI - x(i,4,1))+...
     x(i,4,5).*(VB - x(i,4,1))+...
     x(i,4,4).*(VN - x(i,4,1)).*mg_switch(x(i,4,1)))./CV(4);

%  L5 pyramids: current equations
E = exp(P.H(4,:)).*[0   4   0   0   0   0   0   0] .* m(i,:);
I = exp(P.H(4,:)).*[0   0   0   8   8   0   0   0] .* m(i,:);

E = sum(E) + U(4) + BE;
I = sum(I);

%  L5 pyramids: current equations
f(i,4,2) = ((E' - x(i,4,2)).*KE(i,:)) .* DV(4);
f(i,4,3) = ((I' - x(i,4,3)).*KI(i,:)) .* DV(4);
f(i,4,4) = ((E' - x(i,4,4)).*KN(i,:)) .* DV(4);
f(i,4,5) = ((I' - x(i,4,5)).*KB(i,:)) .* DV(4);

 
% L5 interneurons
%==========================================================================
%  L5 interneurons voltage equations
f(i,5,1) = (GL*(VL - x(i,5,1))+...
     x(i,5,2).*(VE - x(i,5,1))+...
     x(i,5,3).*(VI - x(i,5,1))+...
     x(i,5,5).*(VB - x(i,5,1))+...
     x(i,5,4).*(VN - x(i,5,1)).*mg_switch(x(i,5,1)))./CV(5); 

%  L2/3 interneurons: current equations
E = exp(P.H(5,:)).*[0   0   0   4   0   0   0   0] .* m(i,:);
I = exp(P.H(5,:)).*[0   0   0   0   16  0   0   0] .* m(i,:);

E = sum(E) + BE;
I = sum(I);

%  L5 interneurons: current equations
f(i,5,2) = ((E' - x(i,5,2)).*KE(i,:)) .* DV(5);
f(i,5,3) = ((I' - x(i,5,3)).*KI(i,:)) .* DV(5);
f(i,5,4) = ((E' - x(i,5,4)).*KN(i,:)) .* DV(5);
f(i,5,5) = ((I' - x(i,5,5)).*KB(i,:)) .* DV(5);

 
% L6 thalamic projection pyramids
%==========================================================================
%  L6 pyramids voltage equations
f(i,6,1) = (GL*(VL - x(i,6,1))+...
     x(i,6,2).*(VE - x(i,6,1))+...
     x(i,6,3).*(VI - x(i,6,1))+...
     x(i,6,5).*(VB - x(i,6,1))+...
     x(i,6,4).*(VN - x(i,6,1)).*mg_switch(x(i,6,1)))./CV(6);

%  L6 pyramids: current equations
E = exp(P.H(6,:)).*[0   0   0   2   0   0   0   1/4] .* m(i,:);
I = exp(P.H(6,:)).*[0   0   0   0   8   8   0   0  ] .* m(i,:);

E = sum(E) + BE;
I = sum(I);

%  L6 pyramids: current equations
f(i,6,2) = ((E' - x(i,6,2)).*KE(i,:)) .* DV(6);
f(i,6,3) = ((I' - x(i,6,3)).*KI(i,:)) .* DV(6);
f(i,6,4) = ((E' - x(i,6,4)).*KN(i,:)) .* DV(6);
f(i,6,5) = ((I' - x(i,6,5)).*KB(i,:)) .* DV(6);


% thalamic reticular (inhib)
%==========================================================================
%  thalamic reticular voltage equations
f(i,7,1) = (GL*(VL - x(i,7,1))+...
     x(i,7,2).*(VE - x(i,7,1))+...
     x(i,7,3).*(VI - x(i,7,1))+...
     x(i,7,5).*(VB - x(i,7,1))+...
     x(i,7,4).*(VN - x(i,7,1)).*mg_switch(x(i,7,1)))./CV(7);

%  thalamic reticular: current equations
E = exp(P.H(7,:)).*[0   0   0   0   0   0   0   2] .* m(i,:);
I = exp(P.H(7,:)).*[0   0   0   0   0   0   32  0] .* m(i,:);

E = sum(E) + BE;
I = sum(I);

%  thalamic reticular: current equations
f(i,7,2) = ((E' - x(i,7,2)).*KE(i,:)) .* DV(7);
f(i,7,3) = ((I' - x(i,7,3)).*KI(i,:)) .* DV(7);
f(i,7,4) = ((E' - x(i,7,4)).*KN(i,:)) .* DV(7);
f(i,7,5) = ((I' - x(i,7,5)).*KB(i,:)) .* DV(7);


% thalamic relay (excite)
%==========================================================================
%  thalamic relay voltage equations
f(i,8,1) = (GL*(VL - x(i,8,1))+...
     x(i,8,2).*(VE - x(i,8,1))+...
     x(i,8,3).*(VI - x(i,8,1))+...
     x(i,8,5).*(VB - x(i,8,1))+...
     x(i,8,4).*(VN - x(i,8,1)).*mg_switch(x(i,8,1)))./CV(8);

%  thalamic relay : current equations
E = exp(P.H(8,:)).*[0   0   0   0   0   2   0   0]  .* m(i,:);
I = exp(P.H(8,:)).*[0   0   0   0   0   0   8   32] .* m(i,:);

E = sum(E) + U(1) + BE;
I = sum(I);

%  thalamic relay : current equations
f(i,8,2) = ((E' - x(i,8,2)).*KE(i,:)) .* DV(8);
f(i,8,3) = ((I' - x(i,8,3)).*KI(i,:)) .* DV(8);
f(i,8,4) = ((E' - x(i,8,4)).*KN(i,:)) .* DV(8);
f(i,8,5) = ((I' - x(i,8,5)).*KB(i,:)) .* DV(8);

 
 
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

% Thalamo cortical interactions: 80ms round trip - now in the equations,
% unnecessary here
%--------------------------------------------------------------------------
Tc = zeros(np,np);
%Tc([7 8],[1:6])       = 6; % cortex --> thalamus     = 60 ms
%Tc(1:6,[7 8])         = 2; % thalamus --> cortex     = 20 ms
Tc = -Tc / 1000;
Tc = kron(ones(nk,nk),kron(Tc,eye(ns,ns)));

% Mean intra-population delays, inc. axonal etc. Seem to help oscillation
%--------------------------------------------------------------------------
Dp = ~Ss;                            % states: different sources
Ds = ~Sp & Ss;                       % states: same source different pop.
%Ds = Ds.*(~(Ds & Tc));               % remove t-c and c-t from intrinsic
D  = d(2)*Dp + d(1)*Ds + Tc  ;       %+ Dself;% Complete delay matrix

% Implement: dx(t)/dt = f(x(t - d)) = inv(1 - D.*dfdx)*f(x(t))
%                     = Q*f = Q*J*x(t)
%--------------------------------------------------------------------------
Q  = spm_inv(speye(length(J)) - D.*J);


