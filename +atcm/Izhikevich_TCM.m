function [f,dfdx,Q,D] = Izhikevich_TCM(x,u,P,M)

% x contains hiddent states:
%  x(ns,np,nk): 
%    ns = sources
%    np = cells per source
%    nk = states: 1 = V, 2 = u
%

% states:
x = spm_unvec(x,M.x);
[ns,np,nk] = size(x);

% Input current at t[i]
j = u; u = [];
    
% Loops regions (not really implemented yet)
for i = 1:ns

    I = exp(P.C(1))*j(i);

    % Characteristics of cell types
    RS.a = 0.02; RS.b = 0.2; RS.c = -65; RS.d = 8;    % regular-spiking (RS) cell
    FS.a = 0.1;  FS.b = 0.2; FS.c = -65; FS.d = 2;    % fast-spiking (FS) cell
    TC.a = 0.02; TC.b = 0.25;TC.c = -65; TC.d = 0.05; % thalamo-cortical (TC) cell
    IB.a = 0.02; IB.b = 0.2; IB.c = -55; IB.d = 4;    % intrinsically-bursting (IB) cell

    % Composition of the TC model
    O    = [RS RS FS RS FS IB RS TC];

    % parameters
    a = [O.a].*exp(P.Iz.a);
    b = [O.b].*exp(P.Iz.b);
    c = [O.c];%.*exp(P.Iz.c);
    d = [O.d];%.*exp(P.Iz.d);
    
    a =  a * M.dt;

    V = x(i,:,1);%- 64; 
    %u = b .* V;   %membrane voltage with initial value
    u = x(i,:,2);
    
    tau = M.dt;%0.25;        %tau is the discretization time-step

    for k = 1:8
        if k < 8
            V(k) = V(k) + tau * (0.04 * V(k)^2 + 5 * V(k) + 140 - u(k) );
        else
            V(k) = V(k) + tau * (0.04 * V(k)^2 + 5 * V(k) + 140 - u(k) + I);
        end
    end
    
    %V = V + tau * (0.04 * V.^2 + 5 * V + 140 - u + I); %discretized main equations
    u = u + tau .* a .* (b .* V - u);

    % if there was a spike
    thr = -40;% 30;
    if any( V > thr ) % 30
       isp = find(V>thr);
       mV = V; %VV is the time-series of membrane potentials
       mV(isp) = c(isp);
       %V  = c;
       u(isp)  = u(isp) + d(isp);
       spike_ts = (V > thr); %records a spike
    else
       mV       = V;
       spike_ts = V*0; %records no spike
    end

    
    % compute onnections: if it spikes, it passed on
    H = [...
       -1   0  -1   0   0   1   0   1;
        1  -1  -1   0   0   0   0   0;
        1   1  -1   0   0   0   0   0;
        0   1   0  -1  -1   0   0   0;
        0   0   0   1  -1   0   0   0;
        0   0   0   1  -1  -1   0   1;
        0   0   0   0   0   0  -1   1;
        0   0   0   0   0   1  -1  -1]*4;
    H = H.*exp(P.H(i,:,:));
    
    if any(spike_ts)
        SpikedMat = repmat(spike_ts,[8 1]) .* H;
        SpikedMat = sum(SpikedMat,2);
        %u(find(SpikedMat)) = u(find(SpikedMat))+d(find(SpikedMat)); % make connected neurons sike?
        u = u + SpikedMat';
    end
    
    % revectorise v as output f   
    x(i,:,1) = spm_vec(mV);
    x(i,:,2) = spm_vec(u);
        
end

f = spm_vec(x);

if nargout < 2, return, end

% Only compute Jacobian (gradients) if requested
%==========================================================================
dfdx = spm_cat(spm_diff(M.f,x,u,P,M,1));


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

% Thalamo cortical interactions: ~80ms round trip: 20 ms T->C, 60 ms C->T
%--------------------------------------------------------------------------
%Thalamocortical connections and forward connections from Vp to Vs had
%a mean delay of 3 ms, while corticothalamic connections and backward
%connections from Vs to Vp had a mean delay of 8 m - Neural Dynamics in a Model of the
%Thalamocortical System. I. Layers, Loops and the Emergence of Fast Synchronous Rhythms
% Lumer et al 1997

% CT = 8; %60;
% TC = 3; %20;
% 
% Tc              = zeros(np,np);
% Tc([7 8],[1:6]) = CT  * exp(P.D0(1)); % L6->thal
% Tc([1:6],[7 8]) = TC  * exp(P.D0(2)); % thal->ss
% 
% Tc = -Tc / 1000;
% %Tc = Tc .* ~~(GEa | GIa);
% 
% Tc = kron(ones(nk,nk),kron(Tc,eye(ns,ns)));
% 
% if isfield(P,'ID')
%     % intrisc delays
%     ID = [0 0 0 1 0 1 1 1];
%     ID = [1 .2 .1 1 .2 1 .4 1];
%     ID = [2 1  .1 2 .2 2 .4 2]; % this 
%     %ID = [2 1  1  2  1 2  1 2];
% 
%     ID = -ID.*exp(P.ID)/1000;
%     %ID = -ID/1000;
%     ID = kron(ones(nk,nk),kron(diag(ID),eye(ns,ns)));
% 
%     Tc = Tc + ID;
% end


% if ~isfield(M,'HasThal')
%     M.HasThal = ones(ns,1);
% end
% 
% if all(M.HasThal)
%     Tc = kron(ones(nk,nk),kron(Tc,eye(ns,ns)));
% else
%     Tc = kron(ones(nk,nk),kron(Tc,diag(M.HasThal)));
% end

% Mean intra-population delays, inc. axonal etc. Seem to help oscillation
%--------------------------------------------------------------------------
Dp = ~Ss;                            % states: different sources
Ds = ~Sp & Ss;                       % states: same source different pop.
%Ds = Ds.*(~(Ds & Tc));              % remove t-c and c-t from intrinsic
D  = d(2)*Dp + d(1)*Ds ;%+ Tc  ;       %+ Dself;% Complete delay matrix

% Implement: dx(t)/dt = f(x(t - d)) = inv(1 - D.*dfdx)*f(x(t))
%                     = Q*f = Q*J*x(t)
%--------------------------------------------------------------------------
Q  = spm_inv(speye(length(dfdx)) - D.*dfdx);


end
