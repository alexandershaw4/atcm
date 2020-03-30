function [y,ts,pst] = integratefc(P,M,U,varargin)
% Numerical integration and spectral response of a neural mass model.
% This is the NETWORK version, which handles networks of connected models,
% and with different trial types.
%
% Definitions & methods:
%
% Neural model:
% -------------------------------------------------------------------------
%    dx/dt = f(x,..,P)  // x = states, P = param structure
%
% Numerical integration options:
% -------------------------------------------------------------------------
% (1. ) Straight RK / Euler scheme (if delays are implicit in equations, f):
%
%    y(i+1) = y(i) + dt*f(x,...P)
%
% (2. ) Euler-like scheme with delays:
%
%    dx(t)/dt = f(x(t - d)) = inv(1 - D.*dfdx)*f(x(t))
%                     = Q*f = Q*J*x(t)
% or,
%
%    Q      = spm_inv(speye(length(J)) - D.*J);
%    Q      = (spm_expm(dt*D*J/N) - speye(n,n))*spm_inv(J);
%    y(i+1) = y(i) + Q*dt*f(y(i),..,P)
%
% (3. ) Or a 4th order Runge-Kutta with or without the above delays:
%
%    k1 = f(x          ,...,P)
%    k2 = f(x+0.5*dt*k1,...,P)
%    k3 = f(x+0.5*dt*k2,...,P)
%    k4 = f(x+    dt*k3,...,P)
%
%    y(i+1) = y(i) +     (dt/6) * (k1+2*k2+2*k3+k4)
%    y(i+1) = y(i) + Q * (dt/6) * (k1+2*k2+2*k3+k4)
%
%
% Delays
% -------------------------------------------------------------------------
% Delays between states are dealt with during the integration scheme, using 
% the following system:
%
% Q Implements: dx(t)/dt = f(x(t - d)) = inv(1 - D.*dfdx)*f(x(t))
%                        = Q*f = Q*J*x(t)
%
%
% Observation model:
% -------------------------------------------------------------------------
%    g    :    J'* y + noise              (the weighted signal / LFP timeseries)
%    Pf[i]:    L * ( J[i] * fft(y[i]) )   (over populations, i) 
%
% Integrate & fire scheme:
% -------------------------------------------------------------------------
%    y(i+1) = y(i) + dt*f(x,...,P)
%    S(i)   = 1./(1 + exp(-R.*y(i)')) - 1/2;    [convolution models]
%    F(i)   = sqrt(1 - exp(-(2/pi)*y(i).^2))/2; [conductance models]
%
% Offline, use f2st.m to convert the resultant firing probability time-series 
% to Hz/s series.
%
%
% Dependencies: +atcm tools (which includes:)
%                atcm.Setup.m     - top level set-up script
%                atcm.prepdata.m  - data preparation function
%                atcm.parameters  - prior parameters function
%                atcm.integrate.m - model integration function [this]
%                atcm.tcm.m       - model equations
%                atcm.complete.m  - complete model specification
%                atcm.optim       - tools / functions for optimsation
%                atcm.fun         - misc functions for CSDs etc.
%                atcm.plots       - plotting functions
%
% Also required: SPM12 w/ DCM,  
% AS


% w, initial states, dt | fs is specified & generate sampletimes & inputs
%--------------------------------------------------------------------------
w     = M.Hz;                     % FoI (w)
x     = M.x;                      % model (hidden) states
Kx    = x;                        % pre-fp x, for kernel scheme
dt    = 1/1200;                   % hard-wired 400 hz
Fs    = 1/dt;                     % sampling frequency
tn    = 2;                        % sample window length, in seconds
pst   = 1000*((0:dt:tn-dt)');     % peristim times we'll sample

% unpack simulation pst options
if isfield(M,'sim')
    tn  = M.sim.pst(end);
    pst = M.sim.pst;
    dt  = M.sim.dt;
end

if isfield(M,'intmethod')
    method = M.intmethod;
else
    method = [];
end

if isfield(M,'solvefixed')
    solvefp = M.solvefixed;
else
    solvefp = 1;
end

% Select input type: 0 = constant (DC), 1 = oscilltion, 2 = ERP bump
%--------------------------------------------------------------------------
InputType = 0;

switch InputType
    case 0
        
        % For constant (DC) inputs...
        %------------------------------------------------------------------
        mu    = exp(P.R(1));              % mean amplitude
        drive = ones(length(pst),1)*mu;   % amplitude (constant) over time
        
    case 1
        
        % For oscillatory inputs...
        %------------------------------------------------------------------
        mu    = 1+(P.R(1));                      % mean amplitude
        mf    = 10+(P.R(2));                      % frequency
        drive = mu * sin(2*pi*mf*(pst/1000));  % (sin) oscillation over time
        
    case 2
        
        % For ERP inputs...
        %------------------------------------------------------------------
        delay  = 60 + P.R(1);             % bump
        scale1 = 8  * exp(P.R(2));
        drive  = afit.makef(pst,delay,scale1,16);
        
    case 3
        
        % NOISE
        %------------------------------------------------------------------
        rng default;
        mu    = exp(P.R(1));              % mean amplitude
        hfn   = randn(length(pst),1);
        hfn   = bandpassfilter(hfn,1/dt,[80 300]);
        drive = hfn*mu;   % amplitude (constant) over time
        
end


% expansion (fixed) point: trial & parameter effects are deviations from here
%--------------------------------------------------------------------------
f    = spm_funcheck(M.f); 
if solvefp; x    = atcm.fun.solvefixedpoint(P,M);
else ;      x    = x;%*0;
end

M.x  = x;
v    = spm_vec(x);
NoFX = 0;

if isempty(U)
    U.X  = 1;
    NoFX = 1; % flag no modulatory effects in this model
end

% Integration and spectral response for this trial (c)
%--------------------------------------------------------------------------
for  c = 1:size(U.X,1)
    
    % generate condition-specific parameter structure
    %----------------------------------------------------------------------
    if ~NoFX; Q  = spm_gen_Q(P,U.X(c,:));
    else      Q  = P;
    end
        
    % integration, spectral response, firing etc. (subfunction)
    %----------------------------------------------------------------------
    [y{c} ts{c}] = ...
        dodxdt(pst,f,v,Q,M,dt,w,drive,Kx,U,method,solvefp);

end



end

function [y,ts] = ...
                            dodxdt(t,f,v,P,M,dt,w,drive,Kx,U,method,solvefp)
% Numerical integration, signal weighting and FFT of timeseries with
% spline interpolation and smoothing

% Integration options:
% WithDelays == 0  : Euler without delays
%            == 2  : Euler with delays
%            == 20 : Euler with delays + spm observation function
%            == 5  : Runge-Kutta 45 w/ DCM delays
%            == 45 : Runge-Kutta 45 with delays
%
% Set IntMethod = 'kernels' to override integration and use spm kernels,
% otherwise IntMethod = ' ';
%

if isempty(method)
    method = 2;
end

% Choose either the DCM kernels scheme, or another f(x) update scheme
%--------------------------------------------------------------------------
warning off;
IntMethod   =  ' ';
WithDelays  = method;%2;          % 45 = 4th order Runge-Kutta method with delays
[ns,npp,nk] = size(M.x);

% Prerequisits for integration with delays
if WithDelays == 2 || WithDelays == 5 || WithDelays == 20 || WithDelays == 21 ...
        || WithDelays == 22 || WithDelays == 8 || WithDelays == 24 || WithDelays == 101
    [fx, dfdx,D] = f(M.x,0,P,M);
    
    if any(isnan(fx)) || any(isnan(dfdx(:)))
        fx   = zeros(size(fx));
        dfdx = zeros(size(dfdx));
        D    = zeros(size(D));
    end
    
    OPT.tol = 1e-6*norm((dfdx),'inf');
    
    if OPT.tol == 0
         OPT.tol = 1;
    end
    
    p       = abs(eigs(dfdx,1,'SR',OPT));
    N       = ceil(max(1,dt*p*2));
    n       = spm_length(M.x);
    
    Q       = (spm_expm(dt*D*dfdx/N) - speye(n,n))*spm_inv(dfdx);
    QD      = Q;
    

%     N = 1;
%     %Q = dt*D;  
%     Q  = D*dt;
%     QD = Q;
    
elseif WithDelays == 3
    [fx,~,~,D] = f(M.x,0,P,M);
else
    [fx, dfdx,Q] = f(M.x,0,P,M);
    QD           = Q;
end

try
    N = min([N 4]);
end

% initial firing rate
Curfire = zeros(size(M.x,2),1)';
firings = [];
DoSpecResp = 1;

% % rediscover expansion (fp) point with initial input
if solvefp; M.x    = atcm.fun.solvefixedpoint(P,M);
else ;      M.x    = M.x;%*0;
end
M.pst = t;
v     = spm_vec(M.x);

if WithDelays == 21
    % pre reqs. for int with bilinear jac

    % get Jacobian and its derivatives
    %--------------------------------------------------------------------------
    [dJdx,J]  = spm_diff(f,v,drive(1),P,M,[1 1]);
    [dJdu,J]  = spm_diff(f,v,drive(1),P,M,[1 2]);
    x0        = spm_vec(M.x);
end

        
%fprintf('integrating\n');

% Do an actual numerical integration, if not using kernel approach
%------------------------------------------------------------------
for i   = 1:length(t) % begin time-stepping loop
    
    if ~WithDelays
        % Use a Euler integration scheme
        % y(i+1) = y(i) + dt*f(x,P)
        dxdt   = f(v,drive(i),P,M);
        v      = v + dt*dxdt;
        y(:,i) = v;
        
    elseif WithDelays == 2
        % Karl's Euler-like-with-a-Jacobian-Delay scheme
        % dx = (expm(dt*J) - I)*inv(J)*f(x,u)
        for j = 1:N
            %v = v + Q*f(v,drive(i),P,M,Curfire);
            v = v + Q*f(v,drive(i),P,M);           % CHANGE ME BACK
        end
        % Expansion point - i.e. deviation around fixed point
        y(:,i) = v - spm_vec(M.x);
        
    elseif WithDelays == 20
        % same as above, but with the spm observation (gx) added on
        for j = 1:N
            %v = v + Q*f(v,drive(i),P,M,Curfire);
            v = v + Q*f(v,drive(i),P,M);
        end
        % Expansion about f point
        y (:,i) = v - spm_vec(M.x);
        % Weighted at each integration step!
        yw(:,i) = spm_gx_erp(spm_vec(v),drive(i)',P,M);
        
    elseif WithDelays == 101
        
        dfdx            = spm_diff(M.f,M.x,M.u,P,M,1);
        dfdu            = spm_diff(M.f,M.x,M.u,P,M,2);
        D               = Q;
        dgdx            = spm_diff(M.g,M.x,M.u,P,M,1);
        
        dfdx  = D*dfdx;
        dfdu  = D*dfdu;
        try
            [v0,s] = eig(full(dfdx),'nobalance');
        catch
            v0  = eye(size(dfdx));
            s  = NaN(size(dfdx));
        end
        s     = diag(s);
        
        if max(w) > 1
            s = 1j*imag(s) + real(s) - exp(real(s));
        else
            s = 1j*imag(s) + min(real(s),-1/32);
        end
        
        v       = v + Q*s*dt;
        y (:,i) = v - spm_vec(M.x);
        
        
        
    elseif WithDelays == 21
        % bilinear jacobian
        
        % motion
        %----------------------------------------------------------
        fx    = f(v,drive(i),P,M);
        
        % dx(t)/dt and Jacobian df/dx
        %----------------------------------------------------------
        dx    = spm_vec(v) - x0;
        dfdx  = J;
        for j = 1:length(dJdx)
            dfdx = dfdx + dJdx{j}*dx(j);
        end
        for j = 1:length(dJdu)
            dfdx = dfdx + dJdu{j}*drive(j);
        end
        
        % update dx = (expm(dt*J) - I)*inv(J)*fx
        %----------------------------------------------------------
        v  = spm_unvec(spm_vec(v) + spm_dx(D*dfdx,D*fx,dt),v);
        
        % Treat as expansion about fp
        y(:,i) = v - spm_vec(M.x);;
        
    elseif WithDelays == 22
        
        
        for j = 1:N
            
            % forward Euler:
            %v = v + ( f(v+dt,drive(i),P,M) - f(v,drive(i),P,M) ) / dt;
            
            % backward difference
            %v = v + ( f(v,drive(i),P,M) - f(v+dt,drive(i),P,M) ) / dt;
            
            % central difference
            v = v + Q*( f(v+dt/2,drive(i),P,M) - f(v-dt/2,drive(i),P,M) ) / dt;
            
        end
        
        % Expansion point - i.e. deviation around fixed point
        y(:,i) = v - spm_vec(M.x);
        
        
    elseif WithDelays == 23
        % full jacobian integration
        
        % dx(t)/dt and Jacobian df/dx
        %----------------------------------------------------------------------
        [fx,dfdx,D] = f(v,drive(i),P,M);
        
        % update dx = (expm(dt*J) - I)*inv(J)*fx
        %----------------------------------------------------------------------
        v      = spm_unvec(spm_vec(v) + spm_dx(D*dfdx,D*fx,dt),v);
        
        % output - implement g(x)
        %----------------------------------------------------------------------
        y(:,i) = v - spm_vec(M.x);
        
    elseif WithDelays == 24
        % stochastic equation integration
        
        dfdw = speye(length(v))/sqrt(2);
        
        [fx,dfdx] = f(v,drive(i),P,M);
        v  = spm_unvec(spm_vec(v) + spm_sde_dx(D*dfdx,dfdw,D*fx,dt),v);
        y(:,i) = v - spm_vec(M.x);
        
    elseif WithDelays == 8
        % RK8 !
        
        k_1  = f(v          ,drive(i)               ,P,M);
        k_2  = f(v+dt*(4/27),drive(i)+(dt*4/27)*k_1 ,P,M);
        k_3  = f(v+dt*(2/9) ,drive(i)+  (dt/18)*(k_1+3*k_2),P,M);
        k_4  = f(v+dt*(1/3) ,drive(i)+  (dt/12)*(k_1+3*k_3),P,M);
        k_5  = f(v+dt*(1/2) ,drive(i)+   (dt/8)*(k_1+3*k_4),P,M);
        k_6  = f(v+dt*(2/3) ,drive(i)+  (dt/54)*(13*k_1-27*k_3+42*k_4+8*k_5),P,M);
        k_7  = f(v+dt*(1/6) ,drive(i)+(dt/4320)*(389*k_1-54*k_3+966*k_4-824*k_5+243*k_6),P,M);
        k_8  = f(v+dt       ,drive(i)+  (dt/20)*...
            (-234*k_1+81*k_3-1164*k_4+656*k_5-122*k_6+800*k_7),P,M);
        k_9  = f(v+dt*(5/6) ,drive(i)+ (dt/288)*...
            (-127*k_1+18*k_3-678*k_4+456*k_5-9*k_6+576*k_7+4*k_8),P,M);
        k_10 = f(v+dt       ,drive(i)+ (dt/820)*...
            (1481*k_1-81*k_3+7104*k_4-3376*k_5+72*k_6-5040*k_7-60*k_8+720*k_9),P,M);
        dxdt = dt/840*(41*k_1+27*k_4+272*k_5+27*k_6+216*k_7+216*k_9+41*k_10);
        
        v      = v + Q*dxdt;
        y(:,i) = v - spm_vec(M.x);
        
        
        
    elseif WithDelays == 5
        for j = 1:N
            %4-th order Runge-Kutta method.
            k1       = f(v          ,drive(i)     ,P,M);
            k2       = f(v+0.5*dt*k1,drive(i)+dt/2,P,M);
            k3       = f(v+0.5*dt*k2,drive(i)+dt/2,P,M);
            k4       = f(v+    dt*k3,drive(i)+dt  ,P,M);
            dxdt     = (dt/6)*(k1+2*k2+2*k3+k4);
            v        = v + dxdt;
        end
        y(:,i)   = Q*v;
        
    elseif WithDelays == 45 % RK45 with delayed states
        % 4-th order Runge-Kutta method.
        k1 = f(v          ,drive(i)     ,P,M);
        k2 = f(v+0.5*dt*k1,drive(i)+dt/2,P,M);
        k3 = f(v+0.5*dt*k2,drive(i)+dt/2,P,M);
        k4 = f(v+    dt*k3,drive(i)+dt  ,P,M);
        
        dxdt      = (dt/6)*(k1+2*k2+2*k3+k4);
        v         = v + Q*dxdt;
        y(:,i)    = v - spm_vec(M.x);
        
    end
    
    % firing function at dxdt - assumes conductance model
    %--------------------------------------------------------------
    VR       = -40;
    Vx       = exp(P.S)*32;
    V        = spm_unvec(v,M.x);
    Curfire  = spm_Ncdf_jdw(V(:,:,1),VR,Vx);     % mean firing rate
    S(:,:,i) = Curfire;
    %S = [];
    %fired     = find(squeeze(V(:,:,1)) >= VR);
    %firings   = [firings; [i+0*fired',fired'] ];
    fired=[];
    firings=[];
end

% Reshape to model state space outputs
%--------------------------------------------------------------------------
[ns,npp,nk] = size(M.x);
y(isnan(y)) = 0;
y           = reshape(y,[ns npp nk size(y,2)]);
timeseries  = y;
firing      = S;
nf          = length(w);
spike       = firings;


% Weight - 
J = exp(P.J);
L = exp(P.L);

clear S
for it = 1:length(t)
    for i = 1:length(L)
        S(i,it) = L(i)*sum(sum(J.*squeeze(y(i,:,:,it))));
    end
end


% Compute functional connectivity on weighted -
% for i = 1:ns
%     for j = 1:ns
%         [r(i,j),p(i,j)]=corr( S(i,:)', S(j,:)' );
%     end
% end

band = M.Hz;
S = bandpassfilter(S,1./M.sim.dt,[band(1) band(end)]);
S = abs(hilbert(S));

r=corr(S');
r(isnan(r))=0;

y = kzscorenan(r);
y(isnan(y))=0;
ts=timeseries;
end

function Z=kzscorenan(Z)
s=nanstd(Z(:));
m=nanmean(Z(:));
Z=((Z-m)/s);
end


