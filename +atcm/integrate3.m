function [y,w,s,g,t,pst,layers,noise,firing,QD,Spike] = integrate3(P,M,U,varargin)
% Numerical integration and spectral response of a neural mass model.
% This is the NETWORK version, which handles networks of connected models,
% and with different trial types. This version (atcm.integrate3) uses
% dynamic mode decomposition 
%
% Use as a replacement for spm_csd_mtf.m, for calculating the voltage time-series
% and spectral response from a neural mass.
%
% This routine is more basic than the deault DCM transfer functions because, 
% rather than taking an fft of the system's kernels (eigenvectors of the 
% Jacobian), it uses the  approach of integrating ( dx/dt ) using a Euler 
% or Runge-Kutta scheme over >1s to first generate a the states time series 
% (voltages, currents). Finally, weighting & an fft of  the resultant 
% timeseries provide the spectral response.
%
% This is nice because one can access both the continuous and discrete 
% time responses of the system. Channel noise (exponetial decay functions)
% are computed and returned.
%
% Set DCM.M.IS = @atcm.integrate to use this function with DCM as a 
% replacement for spm_csd_mtf. 
%
% Outputs: y = full model prediction (for comparison to data in DCM.xY.y/csd)
%          w = frequency window (as per DCM.M.Hz)
%          s = the model voltage timeseries for each cell (states)
%          g = the weighted [LFP] timeseries
%          t = the input, over time, to the model
%          pst    = poststim sample times
%          layers = structure with layer specific spectral outputs
%          noise  = structure with (freq space) parameterised noise components 
%          firing = firing rate from each integration step
%
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
% (2.) - THIS VERISON USIING THE OSAKI 1992 METHOD:
%
%    dx(t) = (expm(dfdx*t) - I)*inv(dfdx)*f
%
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
% Observation model: Dynamic mode decomposition
% -------------------------------------------------------------------------
% This version implements dynamic mode decomposition (DMD) and formulates
% the forward projection as a weighted linear combination of discrete
% temporal modes (a Koopman operator). Whereas atcm.integrate weights the 
% states using J, this routine weights the modes using J.
%
% In this version (intgerate3), the DMD routine returns multiple outputs -
% so the system is still 'multiple-output' and quite nonlinear. See
% integrate4.m, where the DMD routine reduces to a single output. This
% converges quicker and prliminary testing suggests might be more accurate.
%
%    yx   :    DMD(y)                     [DMD applied to region state series]
%    Pf[i]:    L * ( J[i] * fft(yx[i]) )  [weighting over modes, i] 
%  c/psd  :    (Pf * Gu) + Gn(auto) + Gs(cross) [augmented w spectral noise]      
%
% Fourier transforms are computed using AfftSmooth.m, which computes cross
% spectra using FFT and the complex-conjugate method. Smoothing is done by
% chunking the series into overlapping windows, computing fft and
% averaging.
%
% Dynamic mode decomposition
% -------------------------------------------------------------------------
%   y           = U*S*V'    
%   y'          = A*U*S*V'
%   U*X'*V*S^-1 = U*AU      = Atilde
%   AtildeW     = WA
%   I           = X'VS^-1*W
%
% Noise model - 3 components (Gu, Gn, Gs)
% -------------------------------------------------------------------------
%   Gu(:,i) = exp(P.a(1,i))*(w.^0);                
%   Gn      = 0;
%   Gs(:,i) = exp(P.c(1,i) - 2)*w.^(-exp(P.c(2,1)));
%
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
%
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

% unpack simulation pst options, if specified
%--------------------------------------------------------------------------
if isfield(M,'sim')
    tn  = M.sim.pst(end);
    pst = M.sim.pst;
    dt  = M.sim.dt;
end

% Select a numerical integration method, if specified
%--------------------------------------------------------------------------
if isfield(M,'intmethod')
    method = M.intmethod;
else
    method = [];
end

% Whether to solve the eqs (f) for a fixed point [default, yes]
%--------------------------------------------------------------------------
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
        mf    = 1+(P.R(2));                      % frequency
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
        hfn   = bandpassfilter(hfn,1/dt,[100 300]);
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
    [y{c},w,s{c},g{c},t{c},layers{c},noise{c},firing{c},QD{c},Spike{c}] = ...
        dodxdt(pst,f,v,Q,M,dt,w,drive,Kx,U,method,solvefp);

end



end

function [y,w,s,g,t,layers,noise,firing,QD,spike] = ...
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
    
    % system delay operator, Q
    Q       = (spm_expm(dt*D*dfdx/N) - speye(n,n))*spm_inv(dfdx);
    QD      = Q;
        
elseif WithDelays == 3
    [fx,~,~,D] = f(M.x,0,P,M);
else
    [fx, dfdx,Q] = f(M.x,0,P,M);
    QD           = Q;
end

try
   % uncomment to put a hard limit on delay integration time
   % N = min([N 4]);
end

% initial firing rate
Curfire = zeros(size(M.x,2),1)';
firings = [];
DoSpecResp = 1;

% % rediscover expansion (fp) point with initial input [no need, done above]
% if solvefp; M.x    = atcm.fun.solvefixedpoint(P,M);
% else ;      M.x    = M.x;%*0;
% end

% Initialise states [v] and period [t] for integration loop
%--------------------------------------------------------------------------
M.pst = t;
v     = spm_vec(M.x);

if WithDelays == 21
    % pre reqs. for int with bilinear jac
    % get Jacobian and its derivatives
    %----------------------------------------------------------------------
    [dJdx,J]  = spm_diff(f,v,drive(1),P,M,[1 1]);
    [dJdu,J]  = spm_diff(f,v,drive(1),P,M,[1 2]);
    x0        = spm_vec(M.x);
end

if isfield(M,'stochastic')
    IsStochastic = M.stochastic;
else
    IsStochastic = 0;
end

% Prespecify matrix & allow the integration to return complex values
y = complex(zeros(ns*npp*nk,length(t)));

switch IntMethod

    case 'kernels'
        
        % Just use the kernels approach?
        %------------------------------------------------------------------
        [K0,K1,K2,H1] = spm_kernels(M,P,length(t),dt);
        y             = H1(:,2:end)';
        S             = [];

        % Transfer functions (FFT of kernel)
        %------------------------------------------------------------------
        S1 = atcm.fun.Afft(K1(:)',1/dt,w);
        DoSpecResp = 0;
        
        % Also compute layer-specific spectra from kernel approach
        %------------------------------------------------------------------
        J  = exp(P.J);
        Ji = find(J);
        for ij = 1:length(Ji)
            P0   = P;
            P0.J = zeros(size(P0.J))-1000;
            P0.J(Ji(ij)) = P.J(Ji(ij));
            [Kj,K1j,Kj,Hj] = spm_kernels(M,P0,length(t),dt);
            Sj(ij,:) = atcm.fun.Afft(K1j(:)',1/dt,w);
        end
        
    case 'ode45'
        
        ode = @(t,v,P,M,f) spm_vec( f(spm_unvec(v,M.x),0,P,M) );
        [~,y]   = ode113(ode,t/1000,spm_vec(v),' ',P,M,f);
        y = y';
        S = [];
        
    otherwise
                
        % Do an actual numerical integration for a discrete epoch, if not using kernel approach
        %------------------------------------------------------------------
        for i   = 1:length(t) % begin time-stepping loop
                        
            if ~WithDelays
                % Use a Euler integration scheme
                % y(i+1) = y(i) + dt*f(x,P)
                dxdt   = f(v,drive(i),P,M);
                v      = v + dt*dxdt;                
                y(:,i) = v;
                                               
            elseif WithDelays == 2 %           [ I RECOMMEND THIS METHOD ]
                
                % Karl's Euler-like-with-a-Jacobian-Delay scheme
                % this just an RK method
                % dx = (expm(dt*J) - I)*inv(J)*f(x,u)
                for j = 1:N
                    %v = v + Q*f(v,drive(i),P,M,Curfire);
                    v = v + Q*f(spm_unvec(v,M.x),drive(i),P,M);           % CHANGE ME BACK   
                    
                    % Ozaki 1992 numerical method:
                    % A bridge between nonlinear time-series models and
                    % nonlinear stochastic dynamical systems: A local 
                    % linearization approach:
                    % "dx(t) = (expm(dfdx*t) - I)*inv(dfdx)*f"
                    %[fx,dfdx] = f(v,drive(i),P,M);
                    %v = v + spm_dx(D*dfdx,D*fx,dt);
                                        
                end   
                % Expansion point - i.e. deviation around fixed point
                if ~IsStochastic
                    y(:,i) = v - spm_vec(M.x);
                else
                    y(:,i) = v - spm_vec(M.x) + rand(size(spm_vec(M.x)));
                end
                
                
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
                
                
            elseif WithDelays == 21
                % Integration with a bilinear jacobian                
                % motion dx(t)/dt and Jacobian df/dx
                %----------------------------------------------------------
                fx    = f(v,drive(i),P,M);
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
                y(:,i) = v - spm_vec(M.x);;  % Treat as expansion about fp          
                                
            elseif WithDelays == 23
                % Full jacobian integration
                % dx(t)/dt and Jacobian df/dx
                %----------------------------------------------------------------------
                [fx,dfdx,D] = f(v,drive(i),P,M);
                
                % update dx = (expm(dt*J) - I)*inv(J)*fx
                %----------------------------------------------------------------------
                v      = spm_unvec(spm_vec(v) + spm_dx(D*dfdx,D*fx,dt),v);
                
                % output - 
                %----------------------------------------------------------------------
                y(:,i) = v - spm_vec(M.x);
                
            elseif WithDelays == 24
                % Stochastic equation integration (SLOW)
                dfdw = speye(length(v))/sqrt(2);
                [fx,dfdx] = f(v,drive(i),P,M);
                v  = spm_unvec(spm_vec(v) + spm_sde_dx(D*dfdx,dfdw,D*fx,dt),v);
                y(:,i) = v - spm_vec(M.x);
                
            elseif WithDelays == 8
                % 8th order Runge Kutta (RK8) - bit ott!
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
            
            
                % firing function at dxdt - assumes conductance model using
                % the JD Williams approximation
                %--------------------------------------------------------------
                VR       = -40;
                Vx       = exp(P.S)*32;
                V        = spm_unvec(v,M.x);
                Curfire  = spm_Ncdf_jdw(V(:,:,1),VR,Vx);     % mean firing rate  
                S(:,:,i) = Curfire;
          %      S = [];
                 %Dfire   = (1./[1 1/2 1/4 1 1/2 1 8 8]);
                 % Dfire   = (1./[1 1   1/4 1 1/2 1 1 1]);   % CHANGE ME BACK!!
                 %Dfire   = Dfire .* exp(P.TV);
                 %Dfire   = (1./[1 1/2 1/4 1 1/2 1 1 1]);
                 %Curfire = Curfire .* Dfire;

                 %y(i) = (-x(i-1) + x(i))./dt;
                 
                %fired     = find(squeeze(V(:,:,1)) >= VR); 
                %firings   = [firings; [i+0*fired',fired'] ];
            fired=[];
            firings=[];
        end
end

warning on;

DoBilinear = 0;
if DoBilinear

    % reduce to a (bi)linear form: operators M0, M1{c}
    %----------------------------------------------------------------------
    [M0,M1,L1,L2] = spm_bireduce(M,P);
    %[M0,M1,L1,L2] = spm_soreduce(M,P);
    M0            = spm_bilinear_condition(M0,length(t),dt);

    % dq/dt = M0*q + u(1)*M1{1}*q + u(2)*M1{2}*q + ....
    %----------------------------------------------------------------------
    M0 = M0(2:end,2:end);    % remove constant
    M1 = M1{1}(2:end,2:end); % remove constant
    qy = M0*y + M1*y;
    %qy = TSNorm(qy,5);
    y  = qy;
    
    % y(i) = L1(i,:)*q + q'*L2{i}*q/2;
    %----------------------------------------------------------------------
    L1 = L1(:,2:end);
    L2 = L2{1}(2:end,2:end);
    for i = 1:length(t)
        dy(:,i)  = L1*qy(:,i) + qy(:,i)'*L2*qy(:,i)/2 ;
        yy(:,i)  =    qy(:,i) + qy(:,i)'*L2*qy(:,i)/2 ;
    end
    %yw          = dy;
    %WithDelays  = 20; % flag to invoke fft(yw)
    y = yy;
end

% Reshape to model state space outputs
%--------------------------------------------------------------------------
[ns,npp,nk] = size(M.x);
y(isnan(y)) = 0;
y(isinf(y)) = 0;
y           = reshape(y,[ns npp nk size(y,2)]);
timeseries  = y;
firing      = S;
nf          = length(w);
spike       = firings;

% Spectral Response Options
%--------------------------------------------------------------------------
%DoSmoothing    = 1; % smooth the fourier series with a 10% Loess window
DoHilbert      = 0; % take the absolute (magnitude) of the hilbert envelope
Bandpassfilter = 0; % band pass filter to [w(1)-1) : w]
DoDCT          = 0; % discrete cosine transform series before fft
IncDCS         = 1; % include a cosine set in the noise model               x       % ON FOR DEXPRO
DoHamming      = 1; % Cosine Hann window the model spectrum      
HamNoise       = 0; % Hann noise components - if data was BPF, exponential 
                    % delay based noise model wont fit without hanning edges
KillTail       = 0;
DoPCA          = 0;
DoSpecResp     = 1; % 1=fft(integrated signal), 2=fft(volterra kernels)

if isfield(M,'IncDCS')
    IncDCS = M.IncDCS;
end


% Compute channel noise components before computing spectral response
%--------------------------------------------------------------------------

% Neuronal innovations: a multiplier on the model signal
%--------------------------------------------------------------------------
for i = 1:size(P.a,2)
    %Gu(:,i) =  exp(P.a(1,i))*(w.^(-exp(P.a(2,i))));
    Gu(:,i) = exp(P.a(1,i))*(w.^0);
end

% Spectrum of channel noise (non-specific): added to spectrum
%--------------------------------------------------------------------------
Gn = P.b(1)*(w.^0)';
%Gn = P.b(1)*w.^(-exp(P.b(2)))';

% Spectrum of channel noise (specific): added to spectrum
%--------------------------------------------------------------------------
for i = 1:size(P.c,2)
    %Gs(:,i) = exp(P.c(1,i) - 2)+w.^(-exp(P.c(2,1)));
    Gs(:,i) = exp(P.c(1,i) )+w.^(-exp(P.c(2,1)));
    %Gs(:,i) = exp(P.c(1,i) - 2)*(w.^0);
    %Gs = Gs*0;
end

%Gs=Gs*0;

if HamNoise
    % Hamming to taper the edges 
    %H  = (1 - cos(2*pi*[1:nf]'/(nf + 1)))/2;
    % optimisation struggles with edge effects
    %----------------------------------------------------------------------
    warning off ;     % dont warn of integrer operands
    H  = kaiser(nf,2.5);
    %Gn = (Gn.*H);
    Gs = (Gs.*H);
    %Gu = (Gu.*H);
    warning on;
    
end

% % Neuronal innovations: a discrete cosine set
%--------------------------------------------------------------------------
if IncDCS
    if isfield(P,'d')
        nd = size(P.d,1);
        X  = spm_dctmtx(nf,nd + 1);
        Mu = exp(X(:,2:end)*P.d);
    else
        Mu = ones(nf,1);
    end
    %if size(Mu,2) == 1, Mu = Mu*ones(1,ns); end       
    %Gu = Gu.*Mu;
    %Gu = exp(P.a(1))*Mu;
    Gu = repmat( (exp(P.a(1))*Mu) , [1 ns]);
end

%[Gu,Gs,Gn] = spm_csd_mtf_gu(P,M.Hz);

% Return noise components for this trial
%--------------------------------------------------------------------------
noise.Gu = Gu;
noise.Gn = Gn;
noise.Gs = Gs;


% Implement the observation model:  y = [ L * fft(J' * y) ] + noise
%--------------------------------------------------------------------------
% Computes the (smoothed) FFT of the weighted (J) signal (g = J*y).
% Adds parameterised channel noise and neuronal innovations in the form of
% a exponential decay functions and a discrete cosine set of order = number
% of populations, respectively.
%
% Use a lengthy burn-in to ensure the system (& optimisation) are in
% steady-state (generating oscillations).
    
% Optional burn in - i.e. take transform of n:end ms instead of 0:end...
burn = findthenearest(300,t); 

% Generate weighted (principal cells) signal - g & fft(g)i
%--------------------------------------------------------------------------
J      = ( exp(P.J(:)) ) ;
Ji     = find(J);
ts     = zeros(ns,length(t));

% This loop generates a weighted 'LFP' time-domain signal for each node in
% the model, with a specific SnR. Note - this is NOT what we'll calculate
% the frequency spectrum from.
for ins = 1:ns
    
    % J-weighted sum of this channel / region
    %----------------------------------------------------------------------
    xseries   = full( squeeze(y(ins,:,:,:)) );
    xseries   = reshape(xseries, [npp*nk,length(t)] );
    ts0 = (sparse(1:8,1,1,size(xseries,1),1)'*xseries)'; % use all mV states
    
    % transform to (smoothed) cosine set if requested
    %----------------------------------------------------------------------
    if DoDCT
        ts0 = dct(ts0);
    end
    
    % add noise to make it look like a real signal
    %----------------------------------------------------------------------
    s2n = mean(ts0(:))*1/16;
    ts0 = ts0 + s2n*randn(size(ts0));
    ts0 = ts0 - mean(ts0);
    
    % apply electrode gain & store this channel / node
    %----------------------------------------------------------------------
    ts(ins,:) = ts(ins,:) + ts0'  * exp(P.L(ins));
    
end

% Select the frequency transfer function - or default fft(dmd(x))
%----------------------------------------------------------------
fmethod = 'dmd'; % svd or dmd

if isfield(M,'fmethod')
    fmethod = M.fmethod;
end

try
switch fmethod
    case 'dmd';
        for ins = 1:ns
            % remove principal [dominant] eigenmode(s) from states series
            x = squeeze( y(ins,:,:,:) );
            x = reshape( x , [npp*nk, length(t)] );
            [u0,s0,v0] = spm_svd((x));
            %ipc = find( (cumsum(diag(full(s0)))./sum(diag(full(s0))) > .8) );
            p1 = u0(:,1)*s0(1,1)*v0(:,1)';
            x = x - p1 ;%- p2;
            %x = full(p1);
            y(ins,:,:,:) = reshape(x,[npp,nk,length(t)]);
        end
        timeseries=y;
end
end

% This is the main loop over regions to calculate the spectrum
%--------------------------------------------------------------------------
for ins = 1:ns
        
    % extract time series of all states from this region
    %----------------------------------------------------------------------
    yx = reshape( squeeze(y(ins,:,:,:)), [npp*nk,length(t)]); 
    yx = yx(1:8,:);        % limit to membrane potentials                    % 1:8 FOR DEXPRO
    %yx = yx([1 2 4 6],:);
    %yx = yx(9:end,:);     % limit to post synaptic currents
    %yx = yx([1 2 4 6],:);    
           
    % parameter controlled suppression of cells: P.Jpop(i)
    if isfield(P,'Jpop') && length(P.Jpop) == 8
        for pij = 1:8
            yx(pij,:) = yx(pij,:)*P.Jpop(pij);
        end
    end
        
    % use dynamic mode decomposition (DMD) to reduce the oscillatory
    % complexity to a simplistic model: a set of discrete spatio-temporal modes
    %----------------------------------------------------------------------
    
    warning off; % avoid rank deficiency warning
        
    %[Eigenvalues, Eigenvectors, ModeAmplitudes, ModeFrequencies, ...
    %    GrowthRates, POD_Mode_Energies] = atcm.fun.dmd(yx, 4, dt);
        
    switch lower(fmethod)
        
%         case 'svd'
%             [Eigenvectors,s,v] = svd(yx);
%             %Eigenvectors = Eigenvectors'*yx;
%             warning on;
% 
%             s  = diag(s);
%             Pf = -(Eigenvectors'*(1./(1j*2*pi*w - s)));
%         
%             Pf = Pf(1:length(Ji),:);
%                 
%             for ij = 1:length(Ji)
%                 for i = 1:length(w)
%                     Pf(ij,i) = sq(Pf(ij,i))*diag(Gu(i,ins))*sq(Pf(ij,i))';
%                 end
%             end
% 
%             for ij = 1:length(Ji)
%                 layers.iweighted (ins,ij,:) = ( Pf(ij,:) * J(ij) ) * exp(P.L(ins));
%             end
%             
%             layers.weighted(ins,:,:) = layers.iweighted(ins,:,:);
            
        case {'dmd' 'instantaneous' 'svd' 'none' 'glm'}
            
            switch fmethod
                case 'dmd'
                %yx = yx([1 2 4 6 8],:);
                [Eigenvalues, Eigenvectors, ModeAmplitudes, ModeFrequencies, ...
                    GrowthRates, POD_Mode_Energies] = atcm.fun.dmd((yx), length(Ji), dt);
                
                %[u,s,v] = svd(yx);
               % 
               % for idmd = 1:4
               %     Eigenvectors(idmd,:) = v(:,idmd)*s(idmd,:)*mean(v,1)';
               % end
                
               %[Eigenvaluesi, Eigenvectorsi, ModeAmplitudesi, ModeFrequenciesi, ...
               %    GrowthRatesi, POD_Mode_Energiesi] = atcm.fun.dmd(real(yx), length(Ji), dt);
                
               % Eigenvectors = Eigenvectors + Eigenvectorsi;
                                
                case 'instantaneous'
                    %[Eigenvectors,s,v] = svd(yx);
                    %Eigenvectors = Eigenvectors'*yx;
                    
                    yx(1,:) = yx(1,:)*0.2;
                    yx(2,:) = yx(2,:)*0.8;
                    yx(3,:) = yx(3,:)*0.1;
                    yx(4,:) = yx(4,:)*0.2;
                    yx(5,:) = yx(5,:)*0.1;
                    yx(6,:) = yx(6,:)*0.2;
                    yx(7,:) = yx(7,:)*0.05;
                    yx(8,:) = yx(8,:)*0.1;
                    %Eigenvectors=yx;
                    [Eigenvectors,s,v] = svd(yx);
                    Eigenvectors = Eigenvectors'*yx;       
                case 'svd'
                    [Eigenvectors,s,v] = svd(cov(real(yx)'));
                    Eigenvectors = Eigenvectors'*yx; 
                case 'none'
                    Eigenvectors = yx;
                case 'glm'
                    % project 8 series onto 4 components, onto 1
                    W = ones(8,4)/(8*4);
                    W = W .* repmat(P.W1(:) ,[1 4]); % 8x1
                    W = W .* repmat(P.W2(:)',[8 1]);
                    
                    Eigenvectors = W'*yx;
                                      
            end
            
            warning on;

            % loop spatio-temporal modes (in this region) and weight them
            %----------------------------------------------------------------------
             for ij = 1:length(Ji)

                y0 = Eigenvectors(ij,:) ;
                Hz = w;

                % Window-smoothed fft 
                %------------------------------------------------------------------
                warning off;
                try
                    
                    switch fmethod
                        case {'none'}
                            % just a smoothed fft of the (contributing)
                            % states
                            [Pf,Hz]  = atcm.fun.AfftSmooth(Eigenvectors(Ji(ij),:),1/dt,w,50);   
                            
                        case {'dmd' 'svd' 'glm'}
                            % just a smoothed fft of the dmd series
                            [Pf,Hz]  = atcm.fun.AfftSmooth(y0(burn:end),1/dt,w,60);          % 60 FOR DEXPRO
                    
                            %[Pf,Hz]  = atcm.fun.Affti(y0(burn:end),1/dt,w); %*
                            %Pf = Pf';
                    
                        case 'instantaneous'
                            % compute the instantaneous frequencies
                            % and amplitudes
                            [ifq,iit] = instfreq(y0(burn:end),1/dt);
                            tburn = t(burn:end)/1000;

                            ufq = unique(round(ifq));
                            rfq = round(ifq);
                            
                            % amplitudes - find time indices
                            for it = 1:length(iit)
                                [~,this(it)] = min(abs(iit(it)-tburn));
                            end
                            
                            % instant amp
                            yh = hilbert(y0);
                            amps = abs(yh(this+burn));
                            
                            % gen gauss bumps
                            spec = w*0;
                            for iq = 1:length(ufq)
                                q  = ufq(iq);
                                %qw = std( ifq(find(ufq(iq)==rfq)) );
                                qw = length( ifq(find(ufq(iq)==rfq)) );
                                
                                amplitude = mean( amps(find(ufq(iq)==rfq)) );
                                
                                pd  = makedist('normal','mu',q,'sigma', 2*(qw));

                                spec = spec + ( abs(amplitude) * pdf(pd,w));

                            end
                            Pf = spec(:)/iq;
                    end
                    
                catch
                    % catch when model goes awol (NaN)
                    Hz = w;
                    Pf = w'*0;
                end
                warning on;

                % Noise shaping - fixed f-scaling
                Pf = (Pf)';
                Pf = full(Pf)';
                Pf = Pf .* Hz';                               % PUT BACK!
                %Pf = abs(Pf);
                
                % Multiply in the parameterised noise terms
                for i = 1:length(Hz)
                    Pf(i,:,:) = sq(Pf(i,:,:))*diag(Gu(i,ins))*sq(Pf(i,:,:))';
                end

                J  = full(J);

                if DoHamming
                    H  = kaiser(nf,2.5);

                    if ~isfield(M,'Hamlower');
                        M.Hamlower = 0;
                    end
                    
                    if M.Hamlower == 0
                        H(1:round(nf/2)) = 1;
                    end

                    Pf = Pf.*H;
                end

                layers.unweighted(ins,ij,:) = ( Pf             )      * exp(P.L(ins));
                layers.weighted  (ins,ij,:) = ( Pf * abs(J(Ji(ij))) ) * exp(P.L(ins));
                layers.iweighted (ins,ij,:) = ( Pf * abs(J(Ji(ij))) ) * exp(P.L(ins));
                layers.DMD(ins,ij,:)        = y0;               % retain DMD series
             end
    end
end

switch lower(fmethod)
    case 'svd'
        clear Pf
end

% Now compute node proper CSDs from sum of weighted cells/modes per region
%----------------------------------------------------------------------
for inx = 1:ns
    for iny = 1:ns
        if ~DoHilbert
            if inx ~= iny
                Pf(:,inx,iny) = sum(layers.iweighted(inx,:,:),2) .* conj( ...
                    sum(layers.iweighted(iny,:,:),2) );
            else
                Pf(:,inx,iny) = sum(layers.iweighted(inx,:,:),2);
            end
        else
            if inx ~= iny
                Pf(:,inx,iny) = max(hilbert(layers.iweighted(inx,:,:)),[],2) .* conj( ...
                    max(hilbert(layers.iweighted(iny,:,:)),[],2) );
            else
                Pf(:,inx,iny) = max(hilbert(layers.iweighted(inx,:,:)),[],2);
            end
        end
    end
end

% Take the absolute (magnitude) of the cross spectra
%----------------------------------------------------------------------
Pf = abs(Pf)/length(w);
%Pf = Pf/length(w);

%    if WithDelays ~= 20
% Incorporate noise components for auto (Gs) and cross (Gn) spectra
%----------------------------------------------------------------------
for i = 1:ns
    %Pf(:,i,i) = Pf(:,i,i) + Gs(:,i);
    for j = 1:ns
        Pf(:,i,j) = Pf(:,i,j) + Gn;
    end
end
        
if DoHamming
    for i = 1:ns
        for j = 1:ns
            H  = (1 - cos(2*pi*[1:nf]'/(nf + 1)))/2;
            H  = kaiser(nf,2.5);
            H(1:round(nf/2)) = 1;
            Pf(:,i,j) = Pf(:,i,j).*H;
        end
    end
end

% kill final portion of frequency window - if data wand bandpass filtered
%--------------------------------------------------------------------------
if KillTail
    % bandpass filteing the real data forces the edges to ~0, but the model
    % can't necessarily do that, so force it here 
    kf   = 70; % tail off from 70
    H    = ones(nf,1); H(end) = 0;
    Ikf  = atcm.fun.findthenearest(w,kf); 
    Wexp = 1./(1:length(Ikf:nf));
    Wexp = smooth(Wexp,.3);
    for i = 1:ns
        for j = 1:ns
            Pf(Ikf:end,i,j) = Pf(Ikf:end,i,j).*Wexp;
        end
    end
end

% returns for this trial - {g}
%--------------------------------------------------------------------------
%y = real(Pf);
y = Pf;
%w = Hz;
s = timeseries;
g = ts;
t = drive;

end

function [x] = sq(x)
% squeeze function for csds
if size(x,3) > 1, x = squeeze(x); else, x = x(:); end

end

function j = adfdx(IS,P,M,order)
% faster computation of the jacobian matrix
% AS
warning off ;

fx    = spm_vec(feval(IS,spm_vec(M.x),0,P,M));
j     = zeros(length(spm_vec(fx)),length(spm_vec(fx)));
for i = 1:length(fx)
                
        % Compute Jacobi: A(j,:) = ( f(x+h) - f(x-h) ) / (2 * h)
        x0     = fx;
        x1     = fx;
        d      = fx(i) * 0.01;
        x0(i)  = x0(i) + d   ;
        x1(i)  = x1(i) - d   ;
                
        f0     = spm_vec(feval(IS,x0,0,P,M));    
        f1     = spm_vec(feval(IS,x1,0,P,M));
        j(i,:) = (f0 - f1) / (2 * d);
        
        if order ==2
        % Alternatively, include curvature
            deriv1 = (f0 - f1) / 2 / d;
            deriv2 = (f0 - 2 * fx + f1) / d ^ 2;
            j(i,:) = deriv1 ./ deriv2;
        end
end

warning on;

j = j';
j(isnan(j))=0;
end

% old delays code:
%             % apply delays to all states
%             FireDel = 1./DV;
%             Df = (FireDel/100)/dt;
%             Df = ceil(Df);
%             for ip = 1:8
%                 if i > Df(ip)
%                     y(ip,i)     = y(ip,i)    - squeeze( y(ip   ,i-(Df(ip)-1) ) );
%                     %y(ip+8,i)  = y(ip+8,i)  + squeeze( y(ip+8 ,i-(Df(ip)-1) ) );
%                     %y(ip+16,i) = y(ip+16,i) + squeeze( y(ip+16,i-(Df(ip)-1) ) );
%                     %y(ip+24,i) = y(ip+24,i) + squeeze( y(ip+24,i-(Df(ip)-1) ) );
%                     %y(ip+32,i) = y(ip+32,i) + squeeze( y(ip+32,i-(Df(ip)-1) ) );
%                 end
%             end
            

%             % manual chunked-fourier
%             j  = sqrt(-1);     % reset j
%             dy = y0(burn:end); % post burn-in
%             FK = 16;           % upsampling
%             fq = resample(w,FK,1);
%             X  = zeros(1,length(fq));
%             
%             NCHUNK = 10;
%             OLW    = round(linspace(1,length(y0),NCHUNK+2)); % overlapping windows
%             
%             for WI  = 1:NCHUNK
%                 dat = y0(OLW(WI):OLW(WI+2));
%                 
%                 for in = 1:length(dat)
%                     for iw = 1:length(fq)
%                         X(iw) = X(iw) + dat(in)*exp(-j*2*pi*(fq(iw)-1)*(in-1)/length(dy));
%                     end
%                 end
%             end
%             
%             for iw = 1:length(w); 
%                 ind(iw) = findthenearest(w(iw),fq);
%             end
%             Pf = abs(X(ind))';
            %Pf = abs(downsample(sX,FK))';