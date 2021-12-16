function [y,w,s,g,t,pst,layers,other] = integrate_1(P,M,U,varargin)
% Numerical integration and spectral response of a neural mass model.
% This is the NETWORK version, which handles networks of connected models,
% and with different trial types. This version (atcm.integrate3) has a number 
% of options for the integration method (Euler, Runge-Kutta 1, 2, 4th) and
% for the transfer function - such as the use of dynamic mode decomposition 
% (DMD), different smothed FFTs and the option to form the output signal
% from a weighted sum of different levels of smoothed time series.
% 
% Usage is designed to be plug-n-play with the way that SPM's Dynamic Causal 
% Modelling (DCM) is structured:
%
% [y,w,s,g,t,pst,layers,noise,firing,QD,Spike] = integrate_1(P,M,U)
%
% Use as a replacement for spm_csd_mtf.m, for calculating the voltage time-series
% and spectral response from a neural mass.
%
% The optimisation problem becomes highly nonlinear so i suggest using AO.m
% optimisation (https://github.com/alexandershaw4/aoptim).
%
% This routine is more basic than the deault DCM transfer functions because, 
% rather than taking an fft of the system's kernels (eigenvectors of the 
% Jacobian), it uses the  approach of integrating ( dx/dt ) using a Euler 
% or Runge-Kutta scheme over >1s to first generate a the states time series 
% (voltages, currents). Finally, weighting & an fft of  the resultant 
% timeseries provide the spectral response.
% This is nice because one can access both the continuous and discrete 
% time responses of the system. Channel noise (exponetial decay functions)
% are computed and returned.
%
% The numerical integration is performed twice: once with no input to the
% system and once with the specified input. The spectrum of each are
% computed and subtracted ensuring that parameter effects are stimulus
% related rather than intrinsic dynamics of the system, but without
% removing induced (ie non-phase locked) effects.
%
% Set DCM.M.IS = @atcm.integrate_1 to use this function with DCM as a 
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
% Definitions, methods and options:
%--------------------------------------
% All options are selected by specifying fields in the input 'M' structure.
%
% Selecting Input Type:
% -------------------------------------
% By default, the input is a constant (D.C).
% Set M.InputType = 0 ... constant
%                 = 1 ... a sine wave oscillation
%                 = 2 ... an ERP (Guassian) bump
%                 = 3 ... high frequency noise
%
% Select Sample Period:
% -------------------------------------
% Default = 2s @ 600Hz.
% Set M.sim.pst = vector of sample times
%     M.sim.dt  = update step (1./fs)
%
% Selecting Numerical Integration Method:
% -------------------------------------
% By default a modified Euler scheme incorporating a delay operator is
% used. 
% To switch from a numerical integration to another option:
% Set M.IntMethod = 'kernels' ... use spm kernel scheme
%       ~~~~      = 'ode45'   ... built in matlab ode solver
%                 = ' '       ... use a numerical method (default)
%                 = 'linearise' ... local linearisation using
%                 eigendecomposition of the jacobian
% 
% If using a numerical method, switch the method type:
% Set M.intmethod = 0  ... Euler no delays
%       ~~~~      = 2  ... Euler with delays
%                 = 21 ... integration with bilinear jacobian
%                 = 23 ... full jacobian integration
%                 = 24 ... stochastic equation integration
%                 = 8  ... 8th-order Runge-Kutta w/ delays
%                 = 45 ... 4th order Runge-Kutta w/ delays
%
% Time series decomposition 
% -------------------------------------
% Set M.decompose = 'ssa'   - singular spectrum analysis algorithm
%                 = 'fourier' - fit a fourier series
%
%
% Selecting Spectral (transfer) Method
% -------------------------------------
% Set M.fmethod = 'dmd' ... use dynamic mode decomposition to identify 
%                           frequency modes in the intgrated timeseries
%               = 'none' ... just use an fft (or smooth fft) of the data
%
% Note - when using 'none', optionally specify whether to take the
% envelope of this spiky spectrum using M.DoEnv (flag 0/1) and how many
% components to include in the envelope (M.ncompe=30) [def=0]
%
% Singular Spectrum Analysis
% -------------------------------------
% The weighted spectral output is decomposed into a frequency basis set 
% using a type of SSA (or EMD or VMD), then projected as a lower dimensional 
% representation composed of the components epxaining the most variance in 
% the full signal (weighted by the data vector in M.y{1}). This helps to 
% 'clean up' noisy spectra in a targetted way.
% 
% Other options:
% -------------------------------------
% M.IncDCS = flag to include a discrete cosine set (a semi-stochastic set
% of neuronal fluctuations in frequency space).
%
% Also required: SPM12 w/ DCM,  
%
% AS

% w, initial states, dt | fs is specified & generate sampletimes & inputs
%--------------------------------------------------------------------------
w     = M.Hz;                     % FoI (w)
x     = M.x;                      % model (hidden) states
Kx    = x;                        % pre-fp x, for kernel scheme
dt    = 1/600;                    % sample rate
Fs    = 1/dt;                     % sampling frequency
tn    = 3;                        % sample window length, in seconds
pst   = 1000*((0:dt:tn-dt)');     % peristim times we'll sample

% unpack simulation pst options, if specified
%--------------------------------------------------------------------------
if isfield(M,'sim')
    pst = M.sim.pst;
    dt  = M.sim.dt;
end

if isfield(P,'dt')
    dt    = (1/1200)*exp(P.dt);
    pst   = 1000*((0:dt:tn-dt)');
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
if isfield(M,'InputType'); InputType=M.InputType;end
switch InputType
    case 0
        % For constant (DC) inputs...
        %------------------------------------------------------------------
        mu    = exp(P.R(1));              % mean amplitude
        drive = ones(length(pst),1)*mu;   % amplitude (constant) over time
        
    case 1
        
        % For oscillatory inputs...
        %------------------------------------------------------------------
        mu    = 1*exp(P.R(1));                      % mean amplitude
        mf    = 1*exp(P.R(2));                      % frequency
        drive = mu * ( (sin(2*pi*mf*(pst/1000))) );%...
                     %   + sin(2*pi*(10*exp(P.R(3)))*(pst/1000)) );
                  
    case 2
        % For ERP inputs...
        %------------------------------------------------------------------
        delay  = 60 * exp(P.R(1));             % bump
        scale1 = 8  * exp(P.R(2));
        drive  = atcm.fun.makef(pst,delay,scale1,16*exp(P.R(3)));
        drive(1)=0;
        
        sust = (max(drive))*.75;
        intcpt = atcm.fun.findthenearest(drive,sust);
        drive(intcpt:end) = sust;%*wave(intcpt:end);
        
    case 3
        % NOISE
        %------------------------------------------------------------------
        rng default;
        mu    = exp(P.R(1));              % mean amplitude
        hfn   = randn(length(pst),1) + (sqrt(-1)*randn(length(pst),1)*1/32);
        drive = hfn*mu;   % amplitude (constant) over time
        drive = .25*drive(1:length(pst));
        
    case 4
        % TWO oscillatory inputs...
        %------------------------------------------------------------------
        mu1   = .001*exp(P.R(1));                % mean amplitude
        mu2   = .001*exp(P.R(2));
        mf1   = 50*exp(P.R(3));                  % frequency
        mf2   = 10*exp(P.R(4));
        
        drive(:,2) = .2*(mu1 * sin(2*pi*mf1*(pst/1000)) ) ;  
        drive(:,1) = .2*(mu2 * sin(2*pi*mf2*(pst/1000)) );
    case 5
        % ifft a parameterised DCT set: neuronal fluctuations
        for i = 1:size(P.a,2)
            Gu(:,i) = exp(P.a(1,i))*(w.^0);        % P.a = constant
        end
        nf = length(w);
        nd = size(P.d,1);
        X  = spm_dctmtx(nf,nd + 1);
        Mu = exp(X(:,2:end)*P.d);
        if size(Mu,2) == 1, Mu = Mu*ones(1,1); end
        Gu = Gu.*Mu;
        Gu = exp(P.a(1))*Mu;

        [Pks,PksAmps]=findpeaks(real(Gu),w);
        for i = 1:length(Pks)
            gx(i,:) = PksAmps(i) * sin(2*pi*Pks(i)*pst);
        end
        drive = sum(gx,1)';
end

if ~isfield(M,'timefreq')
    M.timefreq = 0;
end

if M.timefreq
    if isfield(M,'ons')
        % shift stimulus/drive to onset
        ind = atcm.fun.findthenearest(M.ons,pst);
        newdrive = drive*0;
        newdrive(1:ind)=0;
        newdrive(ind+1:end) = drive(1:length(newdrive(ind+1:end)));
        drive = newdrive;
    end
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
    [y{c},w,s{c},g{c},t{c},layers{c},noise{c},firing{c},QD{c},Spike{c},condel{c},series{c}] = ...
        dodxdt(pst,f,v,Q,M,dt,w,drive,Kx,U,method,solvefp,c);

end

other.noise = noise;
other.firing = firing;
other.QD = QD;
other.Spike = Spike;
other.drive = drive;
other.condel=condel;
other.series = series;
end

function [y,w,s,g,t,layers,noise,firing,QD,spike,condel,series] = ...
                            dodxdt(t,f,v,P,M,dt,w,drive,Kx,U,method,solvefp,ci)
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

if isfield(M,'IntMethod')
    IntMethod = M.IntMethod;
end

% Prerequisits for integration with delays
if WithDelays == 2 || WithDelays == 5 || WithDelays == 20 || WithDelays == 21 ...
        || WithDelays == 22 || WithDelays == 8 || WithDelays == 24 || WithDelays == 101 
    
    % Call the function 
    [fx, dfdx,D] = f(M.x,4,P,M);
    dFdx = dfdx;
    
    % debug instabilities
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
    %Q       = (spm_expm(dt*D*dfdx/N) - speye(n,n))*spm_inv(dfdx);% system delay operator, Q
    %QD      = Q;
    
    dfdx  = dfdx - speye(n,n)*exp(-16);
    Q     = (spm_expm(dt*D*dfdx/N) - speye(n,n))/dfdx;
    
    %Q = sum( ((D.^N).*J)*(Q^N)/factorial(N) ):
    
    QD      = Q;        
elseif WithDelays == 3
    [fx,~,~,D] = f(M.x,0,P,M);
elseif WithDelays == 0
   [fx,dfdx,~,D] = f(M.x,0,P,M);
   QD = D;
else
    [fx, dfdx,Q] = f(M.x,0,P,M);
    QD           = Q;
end

try
   % uncomment to put a hard limit on delay integration time
   %N = min([N 4]);
end

% Set up custom delay operator - parameterised state x state delays
% conferred by population delays (c.f. conduction delay)
%--------------------------------------------------------------------------
if isfield(P,'ID')
    
    if npp == 8 
        del = exp(P.ID).*[2 1/4 1/2 8 1/2 4 2 2]/2.4;
        del = exp(P.ID).*[4 1/4 1 8 1/2 4 2 20]/2.4;
        %ID = [1 1/8 1/4 1 1/2 1 1 8];
        %del = -ID.*exp(P.ID)/1000;
        del = repmat(del,[1 nk]);
        del = 1./del;
        if ns > 1
            if ~isfield(P,'delay')
                del = (spm_vec(repmat(del,[ns 1])))';
            else
            end
        end
        condel=del;
    end
    
else
    del    = 1;
    condel = 1;
end

% initial firing rate
Curfire = zeros(size(M.x,2),1)';
firings = [];
DoSpecResp = 1;

% convert parameterised delay vector to state-by-state matrix
del = sqrt(del)'*sqrt(del);

% limit delays to make computation of derivatives tractable
N = min(N,6);

% Initialise states [v] and period [t] for integration loop
%--------------------------------------------------------------------------
M.pst = t;
v     = spm_vec(M.x);
v0    = v;

if WithDelays == 21
    % pre reqs. for int with bilinear jac-get Jacobian and its derivatives
    %----------------------------------------------------------------------
    [dJdx,J]  = spm_diff(f,v,drive(1),P,M,[1 1]);
    [dJdu,J]  = spm_diff(f,v,drive(1),P,M,[1 2]);
    x0        = spm_vec(M.x);
end

% Frequency steps: dw
dw = 1./(w(2)-w(1));

switch IntMethod

    case 'kernels'
        
        % Just use the kernels approach?
        %------------------------------------------------------------------
        [K0,K1,K2,H1] = spm_kernels(M,P,length(t),dt);
        y = H1(:,1:end-1)';
        S             = [];
        
    case 'ode45'
        % matlab build in ode solver
        ode = @(t,v,P,M,f) spm_vec( f(spm_unvec(v,M.x),0,P,M) );
        [~,y]   = ode113(ode,t/1000,spm_vec(v),' ',P,M,f);
        y = y';
        S = [];
        
    case 'linearise'
        % Here we aim to linearise f, such that:
        % dx = f(x) --> dx = Ax
        % then we use matrix A to predict the timeseries of the system
        
        % initial point
        x0      = spm_vec(M.x);
        % find a fixed point
        if solvefp; xbar    = atcm.fun.solvefixedpoint(P,M);
        else;       xbar    = spm_vec(M.x);
        end
        M.x     = spm_unvec(xbar,M.x);
        
        % compute Jacobian, evaluated at xbar
        dfdx    = spm_diff(M.f,M.x,M.u,P,M,1);
            
        % eigenvectors and values
        [T,D]   = eig(full(real(dfdx)));
        iT      = pinv(T);
        d       = diag(D);
        
        % integrate: x(t) = T*exp(D*t)*iT*x0 
        for i = 1:length(t)
            y(:,i) = T*diag(exp(d*(t(i)./1000)))*iT*xbar;
        end
        S=[];
        
    otherwise
        % Do an actual numerical integration for a discrete epoch, if not using kernel approach
        %------------------------------------------------------------------
        for i   = 1:length(t) % begin time-stepping loop
            if ~WithDelays % (not recommended, more for tinkering)
                % Use a Euler integration scheme
                dxdt   = f(v,drive(i),P,M);
                v      = v + dt*dxdt;  
                y(:,i) = v;
                % With some manually added elays
                Df = real(full(-D)./dt);
                Df = ceil(Df);
                Delays=D*0;
                for ix = 1:size(Df,1)
                    for iy = 1:size(Df,1)
                        if i > abs(Df(ix,iy)) && Df(ix,iy)~=0
                            Delays(ix,iy) = Delays(ix,iy) + y(iy,i-abs(Df(ix,iy)));
                        end
                    end
                end
                if any(Delays(:))
                    Delays = sum(Delays,2)./size(Df,1).^2;
                    y(:,i) = y(:,i) + Delays;
                end
            elseif WithDelays == 2 %           [ I RECOMMEND THIS METHOD ]
                % Karl's Euler-like-with-a-matrix exponential delay operator
                % this just essentially an RK method
                for j = 1:N
                    %v  = v + del'.*(Q*f(spm_unvec(v,M.x),drive(i,:),P,M));
                    %v0 = v0 + del'.*(Q*f(spm_unvec(v0,M.x),0.0001,P,M)); 
                    
                    v  = v + ((del.*Q)*f(spm_unvec(v,M.x),drive(i,:),P,M));
                    v0 = v0 + (del.*Q*f(spm_unvec(v0,M.x),0.0001,P,M)); 
                    
                    %v  = v + (Q*f(spm_unvec(v,M.x),drive(i,:),P,M));
                    %v0 = v0 + (Q*f(spm_unvec(v0,M.x),0.0001,P,M)); 
                end
                y(:,i) = v; %- spm_vec(M.x);  
                y0(:,i) = v0;
                    
            elseif WithDelays == 101
                % logistic equation with harvest as an integration
                r = dt;
                k = 1;
                h = 0;
                for j = 1:N
                    fx = Q*f(spm_unvec(v,M.x),drive(i,:),P,M);            
                    v  = v + r*fx.*( 1 - (fx./k) ) - h;
                end
                y(:,i) = v;
                
            elseif WithDelays == 20
                % same as above, but with the spm observation (gx) added on
                for j = 1:N
                    v = v + Q*f(v,drive(i),P,M);
                end
                % Expansion about f point
                y (:,i) = v - spm_vec(M.x);
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
            Curfire  = spm_Ncdf_jdw(V(:,:,1),VR,Vx);
            S(:,:,i) = Curfire;
            
            % log whether membrane potential crossed threshold
            try
                fired     = find(squeeze(V(:,:,1)) >= VR);
                firings   = [firings; [i+0*fired',fired'] ];
            catch
                % doesn't always work for multi-node models
                fired=[];
                firings=[];
            end
        end
        
end

warning on;

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
try
    yy = spm_unvec(y0,y);
catch
    yy = y;
end

series.States_without = yy;
series.States_with_inp = y;

% EVOKED spectrum: fft(input - no_input)
%[y,s,g,noise,layers] = spectral_response(P,M,y-yy,w,npp,nk,ns,t,nf,timeseries,dt,dfdx,ci,1);

% INDUCED spectrum: ff(input) - fft(no_input)

% system spectral response with specified input
[y1,s,g,noise,layers1] = spectral_response(P,M,y,w,npp,nk,ns,t,nf,timeseries,dt,dfdx,ci,1);

[y0,s,g,noise,layers0] = spectral_response(P,M,yy,w,npp,nk,ns,t,nf,timeseries,dt,dfdx,ci,1);

%y = spm_unvec( (spm_vec(y1) - spm_vec(y0))./spm_vec(y0), y1);
%y = spm_unvec( (spm_vec(y1) - spm_vec(y0)), y1);

%y = atcm.fun.RobustSpectraFit(w,(y1-y0)',1)';

y = y1;

%y(y<0) = -y(y<0);

series.with_inp = y;

% Remove dissipative dynamics leaving only stimulus-related (evoked+induced)  
y = y;%+y0;

% Model the effects of filtering during preprocessing
%--------------------------------------------------------------------------
%if isfield(P,'f')
    %nw = length(w);
    %f     = (1:nw)';
    %f     = exp(real(P.f(1)) + real((P.f(2))*f)/nw);     
%     for i = 1:ns
%         for j = 1:ns
%             %y(:,i,j) = y(:,i,j).*f;
%             y(:,i,j) = atcm.fun.moving_average(y(:,i,j),4);
%         end
%     end
%end


%y = full(exp(P.Ly)*y);
t = drive;

layers = layers1;

%layers = spm_unvec( (spm_vec(layers1)-spm_vec(layers0)), layers1);

end

function [y,s,g,noise,layers]=spectral_response(P,M,y,w,npp,nk,ns,t,nf,timeseries,dt,dfdx,ci,type)

% Spectral Response Options
%--------------------------------------------------------------------------
DoHilbert      = 0; % take the absolute (magnitude) of the hilbert envelope
Bandpassfilter = 0; % band pass filter to [w(1)-1) : w]
DoDCT          = 0; % discrete cosine transform series before fft
IncDCS         = 0; % include semi-stochastic neuronal fluctuations       x       % ON FOR DEXPRO
DoHamming      = 1; %(1)% Cosine Hann window the model spectrum      
HamNoise       = 1; % Hann noise components - if data was BPF, exponential 
                    % delay based noise model wont fit without hanning edges
KillTail       = 0;
DoPCA          = 0;
DoSpecResp     = 1; % 1=fft(integrated signal), 2=fft(volterra kernels)

if isfield(M,'IncDCS')
    IncDCS = M.IncDCS;
end
if isfield(M,'DoHamming')
    DoHamming = M.DoHamming;
end

% Compute channel noise components before computing spectral response
%--------------------------------------------------------------------------
for i = 1:size(P.a,2) % Neuronal innovations: a multiplier on the model signal
    Gu(:,i) = exp(P.a(1,i))*(w.^0);                  % P.a = constant
end

Gn = exp(P.b(1,i) ).*w.^(-exp(P.b(2,1))); % Spectrum of channel noise (non-specific)

for i = 1:size(P.c,2) % Spectrum of channel noise (specific)
    Gs(:,i) = exp(P.c(1,i) )+w.^(-exp(P.c(2,1)));     % P.c = expone
end

if HamNoise
    % Hamming to taper the edges - optimisation struggles with edge effects
    %----------------------------------------------------------------------
    warning off ;        % dont warn of integer operands
    H  = kaiser(nf,2.5);
    warning on;
end

% % Neuronal innovations: a discrete cosine basis set
%--------------------------------------------------------------------------
if IncDCS
    if isfield(P,'d')
        nd = size(P.d,1);
        X  = spm_dctmtx(nf,nd + 1);
        Mu = exp(X(:,2:end)*P.d);
    else
        Mu = ones(nf,1);
    end
    if size(Mu,2) == 1, Mu = Mu*ones(1,ns); end       
    Gu = Gu.*Mu;
    Gu = exp(P.a(1))*Mu;
end

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
burn = atcm.fun.findthenearest(300,t); 

if isfield(M,'burnin')
    burn = atcm.fun.findthenearest(M.burnin,t); 
end

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
    %ts0 = (sparse(1:8,1,1,size(xseries,1),1)'*xseries)'; % use all mV states
    ts0 = spm_vec(J'*xseries);
    
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

% adjustment for when frequency intervals ~= 1 Hz; i.e. dt = dw/dt
%--------------------------------------------------------------------------
dw = 1./(w(2)-w(1));

% This is the main loop over regions to calculate the spectrum
%==========================================================================
for ins = 1:ns
        
    % extract time series of all states from this region
    %----------------------------------------------------------------------
    yx = reshape( squeeze(y(ins,:,:,:)), [npp*nk,length(t)]);
    
    Eigenvectors = yx;
    
    % use dynamic mode decomposition on hidden states if M.dmd = 1
    %----------------------------------------------------------------------
    if isfield(M,'dmd') && M.dmd
        [Eigenvalues, Eigenvectors] = atcm.fun.dmd(yx', 8, dt);
        yx = Eigenvectors*yx;
        Eigenvectors = yx;
        % the user supplied j-vector no longer makes sense, ignore it -
        Ji = 1:8;
        P.J(1:8) = log(1.1);
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
            
            % State timeseries without burnin
            %--------------------------------------------------------------
            this = Eigenvectors(Ji(ij),burn:end);
                        
            % DCT/iDCT components [for filtering only]
            %--------------------------------------------------------------
            [test] = atcm.fun.adct(this)';
            
            pc = test;
            nn = size(test,1);
            
            if all(pc==0)
                test = this;
                pc = test;
                nn = size(test,1);
            end
            
            clear Ppf  Pfm Ppf1 Ppf2 Ppf3
                    
            % Smoothed Fourier transform (atcm.fun.AfftSmooth)
            %--------------------------------------------------------------
            for i = 1:nn; 
                [Ppf(i,:),~,Pfm(i,:,:)] = atcm.fun.AfftSmooth( pc(i,:), dw./dt, w, 30) ;
            end
            
            % Non-smoothed fft version:
            %--------------------------------------------------------------
            %for i = 1:nn; Ppf = atcm.fun.Afft( pc(i,:), dw./dt, w) ;end
            
            % Sliding window smoothing of spectrum
            %--------------------------------------------------------------
            Pf = Ppf;
            Pf = atcm.fun.moving_average(Pf,2);
            
            Pf(isnan(Pf))=0;
            Pf(isinf(Pf))=0;
            
            % return the series / or orthogonal VMD components if selctd
            %--------------------------------------------------------------
            layers.ssa_pc{ins,ij,:,:} = pc;
            layers.pst_burn = t(burn:end);
            
            % Allow a gain specifically for imaginary components
            %--------------------------------------------------------------
            if isfield(P,'iL');
                Pf = real(Pf) + sqrt(-1)*exp(P.iL(ins))*imag(Pf);
            end
            
            Pf  = spm_vec(Pf);
                        
            % Gaussian mixture model - estimate & remove coloured noise
            %-------------------------------------------------------------- 
            tmp = Pf;
            % g = fittype( 'exp(a).*x.^(-exp(b))' );
            c   = atcm.fun.c_oof(w(:),tmp,'exp2');% constrained exponent
            m   = fit(w(:),real(tmp-c),'Gauss4'); % gmm on residual
            Pf  = (c/4) + m(w);

        end
    
        warning on;
        Pf = (Pf)';
        Pf = full(Pf)';
        J  = full(exp(P.J));
        
        if DoHamming
            H = hamming(nf,'periodic');
            Pf = Pf.*H;
        end

        % store the weighted and unweighted population outputs
        %------------------------------------------------------------------
        layers.unweighted(ins,ij,:) = ( Pf             )      * exp(P.L(ins));
        layers.weighted  (ins,ij,:) = ( Pf * abs(J(Ji(ij))) ) * exp(P.L(ins));
        layers.iweighted (ins,ij,:) = ( Pf * abs(J(Ji(ij))) ) * exp(P.L(ins));
    end
end

clear Pf

% Now compute node proper CSDs from sum of weighted cells/modes per region
%--------------------------------------------------------------------------
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
%--------------------------------------------------------------------------
Pf = (Pf)/length(w);

% Hannig window
%--------------------------------------------------------------------------     
if DoHamming
    for i = 1:ns
        for j = 1:ns
            H  = (1 - cos(2*pi*[1:nf]'/(nf + 1)))/2;
            H  = kaiser(nf,2.5);
            Pf(:,i,j) = Pf(:,i,j).*H;
        end
    end
end

%--------------------------------------------------------------------------
if isfield(M,'y')
    for ins = 1:ns
        
        % LFP signal as Weighted sum of population spectra
        %------------------------------------------------------------------
        dat = squeeze(layers.iweighted(ins,:,:))';
        yy  = squeeze(M.y{ci}(:,ins,ins));
        if all( (size(dat,1)==1) && size(dat,2)==length(w) )
            dat = dat';
        end
        Mm = dat';
        b = exp(P.J(Ji));
        %b  = pinv(Mm*Mm')*Mm*yy;
        if any(size(Mm)==length(b))
            Pf(:,ins,ins) = b(:)'*Mm;
            Pf(:,ins,ins) = Pf(:,ins,ins) ;%* exp(P.L(ins));
        else
            Pf(:,ins,ins) = Mm;
        end
                 
        % Electrode gain
        Pf(:,ins,ins) = exp(P.L(ins))*Pf(:,ins,ins);
        
        % return components in separate outputs
        X=[];I=[];RC=[];
        b=[];pci=[];
        c = pci;X = b;
        meanpower={X c};
            
        % Add noise to (in frequency space) to this LFP channel spectrum
        %------------------------------------------------------------------
        addnoise=1;
        if addnoise
            % Multiply in the {semi-stochastic} neuronal fluctuations
            % note this would be a convolution of signal and noise in the time domain
            for i = 1:length(Hz)
                Pf(i,ins,ins) = (sq(real(Pf(i,ins,ins)))*real(diag(Gu(i,ins)))*sq(real(Pf(i,ins,ins)))) + sqrt(-1)*imag(Pf(i,ins,ins));
            end
        end
        
    end                                            % end loop over sources

    % Recompute cross spectrum of regions (CSDs)s
    %------------------------------------------------------------------
    for inx = 1:ns
        for iny = 1:ns
            if inx ~= iny
                Pf(:,inx,iny) = squeeze( Pf(:,inx,inx) ) .* conj( Pf(:,iny,iny) ) ;
            end
        end
    end
    
    % Re-incorporate other noise components for auto (Gs) and cross (Gn) spectra
    %----------------------------------------------------------------------
    %if ns > 1
    addnoise=1;
    if addnoise        
        for i = 1:ns
            for j = 1:ns
                % Autospectral noise / innovations [P.b]
                Pf(:,i,j) = (Pf(:,i,j).*Gs(:,i)) + Gn(:);
                if j ~= i
                    % Cross spectral noise / innovations [P.c]
                    Pf(:,j,i) = Pf(:,j,i) .* Gs(:,i);
                    Pf(:,i,j) = Pf(:,j,i);
                end
            end
        end
     end
    
end

% Model the effects of filtering during preprocessing
%--------------------------------------------------------------------------
if isfield(P,'f')
    nw = length(w);
    f     = (1:nw)';
    f     = exp(real(P.f(1)) + real((P.f(2))*f)/nw);     
    for i = 1:ns
        for j = 1:ns
            Pf(:,i,j) = Pf(:,i,j).*f;
        end
    end
end

% returns for this trial - {g}
%--------------------------------------------------------------------------
y = Pf;
s = timeseries;
g = ts;

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
