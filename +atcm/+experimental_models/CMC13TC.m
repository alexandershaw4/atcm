function [f,J,Q] = CMC13TC(x,u,P,M)
% state equations for a neural mass model (canonical microcircuit)
% FORMAT [f,J,D] = spm_fx_cmc(x,u,P,M)
% FORMAT [f,J]   = spm_fx_cmc(x,u,P,M)
% FORMAT [f]     = spm_fx_cmc(x,u,P,M)
% x      - state vector
%
%   x(:,1) - voltage     (spiny stellate cells)
%   x(:,2) - conductance (spiny stellate cells)
%   x(:,3) - voltage     (superficial pyramidal cells)
%   x(:,4) - conductance (superficial pyramidal cells)
%   x(:,5) - current     (inhibitory interneurons)
%   x(:,6) - conductance (inhibitory interneurons)
%
%   x(:,7)  - voltage     (deep L5 pyramidal cells)
%   x(:,8)  - conductance (deep L5 pyramidal cells)
%   x(:,9)  - voltage     (deep interneurons)
%   x(:,10) - conductance (deep interneurons)
%   x(:,11) - voltage     (deep L6 Pyramidal cells)
%   x(:,12) - conductance (deep L6 Pyramidal cells)

%   x(:,13) - voltage     (thalamic reticular cells)
%   x(:,14) - conductance (thalamic reticular cells)
%   x(:,15) - voltage     (thalamic relay cells)
%   x(:,16) - conductance (thalamic relay cells)
%
%
% f        - dx(t)/dt  = f(x(t))
% J        - df(t)/dx(t)
% D        - delay operator dx(t)/dt = f(x(t - d))
%                                    = D(d)*f(x(t))
%
% Prior fixed parameter scaling [Defaults]
%
% E  = (forward, backward, lateral) extrinsic rates 
% G  = intrinsic rates
% D  = propagation delays (intrinsic, extrinsic)
% H  = overall synaptic kinetics
% T  = synaptic time constants
% R  = slope of sigmoid activation function
%
%__________________________________________________________________________
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%___________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_fx_cmc.m 4348 2011-06-10 20:50:23Z karl $
 
 
% get dimensions and configure state variables
%--------------------------------------------------------------------------
M.u   = u;                          % place inputs in M
x     = spm_unvec(x,M.x);           % neuronal states
[n m] = size(x);                    % number of sources and states  
 
 
% [default] fixed parameters
%--------------------------------------------------------------------------
E  = [1 1/2 1 1/2]*200;                % extrinsic (forward and backward)  
%G  = [4 4 4 4 4 2 4 4 2 1 4 4 4]*200; % intrinsic connections 
D  = [1 16];                           % delays (intrinsic, extrinsic)
T  = [2 2 10   20 10 20   10 20];             
R  = 2/3;                              % slope of sigmoid activation function

% Alex: additional G's for new thalamic connections
%G = [G [2 2 2 4 1 1]*200];
%T = [T 60 20];
 

G = [...
... ss  sp  si  dp  di  tp  rt  rl   
    4   0   4   0   0   4   0   4;  
    4   4   4   0   0   0   0   0;  
    4   4   4   0   0   0   0   0;
    0   4   0   4   4   0   0   0;
    0   0   0   4   4   0   0   0;
    0   0   0   4   4   4   0   4;
    0   0   0   0   0   0   4   4;
    0   0   0   0   0   4   4   4 ];

G = G.*exp(P.G);

G0 = [...
... ss  sp  si  dp  di  tp  rt  rl   
   -1   0  -1   0   0   1   0   1;  
    1  -1  -1   0   0   0   0   0;  
    1   1  -1   0   0   0   0   0;
    0   1   0  -1  -1   0   0   0;
    0   0   0   1  -1   0   0   0;
    0   0   0   1  -1  -1   0   1;
    0   0   0   0   0   0  -1   1;
    0   0   0   0   0   1   1  -1 ];

 
% Extrinsic connections
%--------------------------------------------------------------------------
% ss = spiny stellate
% sp = superficial pyramidal
% dp = deep pyramidal
% ii = inhibitory interneurons
%--------------------------------------------------------------------------
A{1} = exp(P.A{1})*E(1);           % forward  connections (sp -> ss)
A{2} = exp(P.A{2})*E(2);           % forward  connections (sp -> dp)
A{3} = exp(P.A{3})*E(3);           % backward connections (dp -> sp)
A{4} = exp(P.A{4})*E(4);           % backward connections (dp -> ii)
C    = exp(P.C);
 
% pre-synaptic inputs: s(V)
%--------------------------------------------------------------------------
R    = R.*exp(P.S);
S    = 1./(1 + exp(-R*x)) - 1/2;
 
% exogenous input
%--------------------------------------------------------------------------
U    = C*u(:);
 

% time constants and intrinsic connections
%==========================================================================
T    = ones(n,1)*T/1000;

T = T.*exp(P.T);

%G    = ones(n,1)*G;
 
% free parameters on time constants and intrinsic connections
%--------------------------------------------------------------------------
% G(:,1)  ss -> ss
% G(:,2)  sp -> ss
% G(:,3)  ii -> ss
% G(:,4)  ii -> ii
% G(:,5)  ss -> ii
% G(:,6)  dp -> ii
% G(:,7)  sp -> sp
% G(:,8)  ss -> sp
% G(:,9)  ii -> dp
% G(:,10) dp -> dp
% G(:,11) ii -> sp ADDED BY KRISH
% G(:,12) sp -> ii ADDED BY KRISH
% G(:,13) sp -> dp ADDED BY KRISH
%
% G(:,14) dp -> rt [added by alex]
% G(:,15) rt -> rl [added by alex]
% G(:,16) rl -> rt [added by alex]
% G(:,17) rl -> ss [added by alex]
% G(:,18) rt -> rt
% G(:,19) rl -> rl
%
%--------------------------------------------------------------------------

% ModT=P.T;
%     for i=1:size(P.T,1);
%         for j=2:size(P.T,2);
%             if ModT(i,j)==-100.
%                 ModT(i,j)=ModT(i,1);
%             end
%         end
%     end
% 
% for i = 1:size(P.T,2)
%     T(:,i) = T(:,i).*exp(ModT(:,i));
% end
% for i = 1:size(P.G,2)
%     G(:,i) = G(:,i).*exp(P.G(:,i));
% end
 
%f = x;%*0;

% Motion of states: f(x)
%--------------------------------------------------------------------------
 
% Conductance
%==========================================================================
%
% S(:,1) --- Stellates ss
% S(:,3) --- Superficial Pyramidals sp
% S(:,5) --- Inibitory Interneurons ii
% S(:,7) --- Deep pyramidals dp
% S(:,9) --- Thalamic Reticular
% S(:,11) -- Thalamic Relay
%
%

% Prepare cortico-cortical and ex/endogenous currents
uc{1} = + A{1}*S(:,3);
uc{2} = - A{3}*S(:,7);
uc{3} = - A{4}*S(:,7);
uc{4} = + A{2}*S(:,3);
uc{5} = 0;
uc{6} = 0;
uc{7} = 0;
uc{8} = U;

% Current state indices for the 8 pops:
ic = 2:2:16;

% Equations from matrix form
for i = 1:8
    
    u    = 0;
    from = find( G(i,:) );
    
    for j = 1:length(from)
        u = u + ( G0(i,from(j)) * G(i,from(j)).*S(:,from(j)) );
    end
    
    f(:,ic(i)) = (u - 2*x(:,ic(i)) - x(:,i)./T(:,i))./T(:,i);
    
end

% Voltage
for i = 1:length(ic)
    f(:,ic(i)-1) = f(:,ic(i));
end


% % Granular layer (excitatory interneurons): spiny stellate: Hidden causes
% %--------------------------------------------------------------------------
% u      =   A{1}*S(:,3) + U;
% u      = - G(:,1).*S(:,1) - G(:,3).*S(:,5) - G(:,2).*S(:,3) + u;
% u      =   G(:,17).*S(:,11) + u;
% f(:,2) = (u - 2*x(:,2) - x(:,1)./T(:,1))./T(:,1);
%  
% % Supra-granular layer (superficial pyramidal cells): Hidden causes - error
% %--------------------------------------------------------------------------
% u      = - A{3}*S(:,7);
% u      =   G(:,8).*S(:,1) - G(:,7).*S(:,3) + u;
% u      = - G(:,11).*S(:,5) + u; % KRISH ADDED ii -> sp modulation.
% f(:,4) = (u - 2*x(:,4) - x(:,3)./T(:,2))./T(:,2);
%  
% % Supra-granular layer (inhibitory interneurons): Hidden states - error
% %--------------------------------------------------------------------------
% u      = - A{4}*S(:,7);
% u      =   G(:,5).*S(:,1) + G(:,6).*S(:,7) - G(:,4).*S(:,5) + u;
% u      =  G(:,12).*S(:,3) + u; % KRISH ADDED sp -> ii modulation.
% f(:,6) = (u - 2*x(:,6) - x(:,5)./T(:,3))./T(:,3);
%  
% % Infra-granular layer (deep pyramidal cells): Hidden states
% %--------------------------------------------------------------------------
% u      =   A{2}*S(:,3);
% u      = - G(:,10).*S(:,7) - G(:,9).*S(:,5) + u;
% u      =  G(:,13).*S(:,3) + u; % KRISH ADDED sp -> dp modulation.
% f(:,8) = (u - 2*x(:,8) - x(:,7)./T(:,4))./T(:,4);
%  
% % Thalamic RETICULAR cells
% %-------------------------------------------------------------------
% u       = G(:,14).*S(:,7);      % dp -> rt
% u       = G(:,16).*S(:,11) + u; % rl -> rt
% u       = -G(:,18).*S(:,9) + u; % rt -> rt
% f(:,10) = (u - 2*x(:,10) - x(:,9)./T(:,5))./T(:,5);
% 
% % Thalamic RELAY cells
% %-------------------------------------------------------------------
% u       = -G(:,15).*S(:,9);      % rt -> rl
% u       = -G(:,19).*S(:,11) + u; % rl -> rl
% f(:,12) = (u - 2*x(:,12) - x(:,11)./T(:,6))./T(:,6);

% % Voltage
% %==========================================================================
% f(:,1)  = x(:,2);
% f(:,3)  = x(:,4);
% f(:,5)  = x(:,6);
% f(:,7)  = x(:,8);
% f(:,9)  = x(:,10);
% f(:,11) = x(:,12);
f       = spm_vec(f);

%G
%S
%f

if nargout == 1; return, end
 
 
% delays
%==========================================================================
% Delay differential equations can be integrated efficiently (but
% approximately) by absorbing the delay operator into the Jacobian
%
%    dx(t)/dt     = f(x(t - d))
%                 = Q(d)f(x(t))
%
%    J(d)         = Q(d)df/dx
%--------------------------------------------------------------------------
De = exp(P.D);
Di = diag(diag(De));
De = De - Di;
De = De*D(2)/1000;
Di = Di*D(1)/1000;
De = kron(ones(m,m),De);
Di = kron(ones(m,m) - speye(m,m),Di);
D  = Di + De;

% D(1:8,9:12) = 20.*exp(P.D(2)); % thal -> cort
% D(9:12,1:8) = 60.*exp(P.D(3)); % cort -> thal

 
% Implement: dx(t)/dt = f(x(t - d)) = inv(1 + D.*dfdx)*f(x(t))
%                     = Q*f = Q*J*x(t)
%--------------------------------------------------------------------------
try 
    [Q,J] = spm_dcm_delay(M,P,D);
catch
    [Q,J] = spm_dcm_delay(P,M);
end
 
 
return
 
% notes and alpha function (kernels)
%==========================================================================
% x   = t*exp(k*t)
% x'  = exp(k*t) + k*t*exp(k*t)
%     = exp(k*t) + k*x
% x'' = 2*k*exp(k*t) + k^2*t*exp(k*t)
%     = 2*k*(x' - k*x) + k^2*x
%     = 2*k*x' - k^2*x
