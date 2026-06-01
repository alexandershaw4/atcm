function [f,J,Q] = spm_fx_venst(x,u,P,M)
% state equations for  5 pop ventral striatum
%
% FORMAT [f,J,Q] = spm_fx_tcm(x,u,P,M)
%
% x - states and covariances
%
% x(i,j,k)        - k-th state of j-th population of i-th source
%                   i.e., running over sources, pop. and states
%
%   population: 1 - excitatory spiny stellate cells (input cells)
%               2 - superficial pyramidal cells     (forward output cells)
%               3 - inhibitory interneurons         (intrisic interneurons)
%               4 - deep pyramidal cells            (backward output cells)
%               5 - deep interneurons
%               6 - Thalamic projection pyramidal neurons (L6)
%
%        state: 1 V  - voltage
%               2 gE - conductance (excitatory)
%               3 gI - conductance (inhibitory)
%
%
% Note: This is the most vanilla implementation of TCM, which seems to work
% well for CSD models using the CSD tranfer function (i.e. spm_csd_mtf.m),
% but less with numerical integration methods, such as the euler and RK methods 
% implemented in atcm.integrate_1. this version of TCM is just an extended
% version of spm_fx_cmm_NMDA.m but with the addition of the extra
% populations and a new (parameterised in P.D0) delay operator.
%
% {AS2021}
%--------------------------------------------------------------------------
% refs:
%
% Marreiros et al (2008) Population dynamics under the Laplace assumption
%
% See also:
%
% Friston KJ.
% The labile brain. I. Neuronal transients and nonlinear coupling. Philos
% Trans R Soc Lond B Biol Sci. 2000 Feb 29;355(1394):215-36. 
% 
% McCormick DA, Connors BW, Lighthall JW, Prince DA.
% Comparative electrophysiology of pyramidal and sparsely spiny stellate
% neurons of the neocortex. J Neurophysiol. 1985 Oct;54(4):782-806.
% 
% Brunel N, Wang XJ.
% What determines the frequency of fast network oscillations with irregular
% neural discharges? I. Synaptic dynamics and excitation-inhibition
% balance. J Neurophysiol. 2003 Jul;90(1):415-30.
% 
% Brunel N, Wang XJ.
% Effects of neuromodulation in a cortical network model of object working
% memory dominated by recurrent inhibition. J Comput Neurosci. 2001
% Jul-Aug;11(1):63-85.
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_fx_cmm_NMDA.m 5741 2013-11-13 12:10:48Z guillaume $
 
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
A{1}  = exp(P.A{1});                      % forward
A{2}  = exp(P.A{2});                      % backward
AN{1} = exp(P.AN{1});                      % forward
AN{2} = exp(P.AN{2});                      % backward
C     = exp(P.C);                         % subcortical
 

% detect and reduce the strength of reciprocal (lateral) connections
%--------------------------------------------------------------------------
for i = 1:length(A)
    L    = (A{i} > exp(-8)) & (A{i}' > exp(-8));
    A{i} = A{i}./(1 + 8*L);
end

            
% intrinsic connection strengths
%==========================================================================

% condition specific effects: Inhibition of SPGP's
%--------------------------------------------------------------------------
G    = full(P.H);
if any(P.G)
    G(2,2,:) = squeeze(G(2,2,:)) + P.G;
end
  G    = exp(G);



% connectivity switches
%==========================================================================
% 1 - excitatory spiny stellate cells (granular input cells)
% 2 - superficial pyramidal cells     (forward  output cells)
% 3 - inhibitory interneurons         (intrisic interneuons)
% 4 - deep pyramidal cells            (backward output cells)


% extrinsic connections (F B) - from superficial and deep pyramidal cells
%--------------------------------------------------------------------------
SA   = [1   0 ; % spn 1
        0   0 ; % sst
        0   0 ; % ii  1
        0   1 ; % spn 2
        0   0]; % ii  2
% extrinsic NMDA-mediated connections (F B) - from superficial and deep pyramidal cells
%--------------------------------------------------------------------------    
    
SNMDA   = [1   0 ;
           0   0 ;
           0   0 ;
           0   1 ;
           0   0];
    

% intrinsic connections (np x np) - excitatory
%--------------------------------------------------------------------------
GE   = [ 0    1    0    0    0; % spn 1
         1    0    0    1    0; % sst 
         1    0    0    1    0; % ii
         0    1    0    0    0; % spn 2
         0    0    0    0    0];% ii 2
     
% intrinsic connections (np x np) - inhibitory
%--------------------------------------------------------------------------
GI   = [ 0    0    1    0    0; % spn 1
         0    0    0    0    0; % sst 
         0    0    0    0    1; % ii
         0    0    1    0    0; % spn 2
         0    0    1    0    1];% ii 2
     
% rate constants (ns x np) (excitatory 4ms, inhibitory 16ms)
%--------------------------------------------------------------------------
KE    = exp(-P.T(:,1))*1000/4;                       % excitatory rate constants (AMPA)
KI    = exp(-P.T(:,2))*1000/16;                      % inhibitory rate constants (GABAa)
KNMDA = exp(-(P.T(:,3) ))*1000/100;                  % excitatory rate constants (NMDA)


% Voltages
%--------------------------------------------------------------------------
VL   = -70;                               % reversal  potential leak (K)
VE   =  60;                               % reversal  potential excite (Na)
VI   = -90;                               % reversal  potential inhib (Cl)
VR   = -40;                               % threshold potential
VN   =  10;                               % reversal Ca(NMDA)   

CV   = exp(P.CV).*[128 128 256 128 256]/1000;  % membrane capacitance
GL   = 1;                                 % leak conductance
 
% mean-field effects:
%==========================================================================

% neural-mass approximation to covariance of states: trial specific
%----------------------------------------------------------------------
Vx   = exp(P.S)*32;

% mean population firing and afferent extrinsic input
%--------------------------------------------------------------------------
  
m       = spm_Ncdf_jdw(x(:,:,1),VR,Vx);     % mean firing rate  
a(:,1)  = A{1}*m(:,2);                      % forward afference  AMPA
a(:,2)  = A{2}*m(:,4);                      % backward afference AMPA 
an(:,1) = AN{1}*m(:,2);                     % forward afference  NMDA
an(:,2) = AN{2}*m(:,4);                     % backward afference NMDA

% Averge background activity and exogenous input
%==========================================================================
BE     = exp(P.E)*0.8;

% input
%--------------------------------------------------------------------------
if isfield(M,'u')
    
    % endogenous input
    %----------------------------------------------------------------------
    U = u(:);
    
else
    
    % exogenous input
    %----------------------------------------------------------------------
    U = C*u(:);
    
end

% flow over every (ns x np) subpopulation
%==========================================================================
f     = x;

for i = 1;%:ns
   
        ip = 1:5;
    
        % intrinsic coupling
        %------------------------------------------------------------------
        E      = (G(:,:,1).*GE)*m(2,ip)';
        ENMDA  = (G(:,:,1).*GE)*m(2,ip)';
        I      = (G(:,:,1).*GI)*m(2,ip)';
        
        
        % extrinsic coupling (excitatory only) and background activity
        %------------------------------------------------------------------
        E     = (E +  BE + SA*a(2,:)')*2;
        ENMDA = (ENMDA + BE  + SNMDA*an(2,:)')*2;
      

        % and exogenous input(U)
        %------------------------------------------------------------------
        E(1) = E(1)  + U(i);
        ENMDA(1) = ENMDA(1) + U(i);
        
        % Voltage
        %==================================================================
          f(i,ip,1) =    (GL*(VL - x(i,ip,1))+...
                         x(i,ip,2).*(VE - x(i,ip,1))+...
                         x(i,ip,3).*(VI - x(i,ip,1))+...
                         x(i,ip,4).*(VN - x(i,ip,1)).*mg_switch(x(i,ip,1)))./CV;
        
        % Conductance
        %==================================================================
        f(i,ip,2) = (E' - x(i,ip,2)).*KE(i,:);
        f(i,ip,3) = (I' - x(i,ip,3)).*KI(i,:);
        f(i,ip,4) = (ENMDA' - x(i,ip,4))*KNMDA(i,:) ;
end
          
           
% vectorise equations of motion
%==========================================================================
f = spm_vec(f);
 
if nargout < 2, return, end

% Jacobian
%==========================================================================
J = spm_cat(spm_diff(M.f,x,u,P,M,1));

if nargout < 3, return, end

% Delays
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

d  = -D.*exp(P.D)/1000;
Sp = kron(ones(nk,nk),kron( eye(np,np),eye(ns,ns)));  % states: same pop.
Ss = kron(ones(nk,nk),kron(ones(np,np),eye(ns,ns)));  % states: same source

Dp = ~Ss;                            % states: different sources
Ds = ~Sp & Ss;                       % states: same source different pop.
D  = d(2)*Dp + d(1)*Ds;         % add them including thalamo-cortical 


% Implement: dx(t)/dt = f(x(t - d)) = inv(1 - D.*dfdx)*f(x(t))
%                     = Q*f = Q*J*x(t)
%--------------------------------------------------------------------------
Q  = spm_inv(speye(length(J)) - D.*J);


