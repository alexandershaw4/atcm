function [y,w,s,g,t,pst,layers,other] = integrate3(P,M,U,varargin)
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
% [y,w,s,g,t,pst,layers,noise,firing,QD,Spike] = integrate3(P,M,U)
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
% This is nice because one can access both the continuous and discrete 
% time responses of the system. Channel noise (exponetial decay functions)
% are computed and returned.
%
% Set DCM.M.IS = @atcm.integrate3 to use this function with DCM as a 
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
% Selecting Spectral (transfer) Method
% -------------------------------------
% Set M.fmethod = 'dmd' ... use dynamic mode decomposition to identify 
%                           frequency modes in the intgrated timeseries
%               = 'instantaneous' ... estimate instantaneous frequency and
%                            amplitude from the phase of the int'd series
%               = 'svd'  ... use an svd-based PCA to idenfiy frequencies in
%                            integrated timseries
%               = 'none' ... just use an fft (or smooth fft) of the data
%               = 'timefreq' ... do a time-frequency decomp and average it
%
% Note - when using 'none', optionally specify whether to take the
% envelope of this spiky spectrum using M.DoEnv (flag 0/1) and how many
% components to include in the envelope (M.ncompe=30).
%
% Other options:
% -------------------------------------
% M.IncDCS = flag to include a discrete cosine set (a semi-stochastic set
% of neuronal fluctuations in frequency space).
%
% Including the data (DCM.xY.y) in M (i.e. DCM.M.y = DCM.xY.y) will invoke
% a linear model, formed by a combination of the raw fourier spectrum and
% the corresponding smoothed (envelope) version. This 
% is a way of 'optimising' the smoothing function.
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
        mu    = 1*(P.R(1));                      % mean amplitude
        mf    = 10*exp(P.R(2));                      % frequency
        drive = mu * ( (sin(2*pi*mf*(pst/1000))) );%...
                     %   + sin(2*pi*(10*exp(P.R(3)))*(pst/1000)) );
                    
    case 2
        % For ERP inputs...
        %------------------------------------------------------------------
        delay  = 60 + P.R(1);             % bump
        scale1 = 8  * exp(P.R(2));
        drive  = atcm.fun.makef(pst,delay,scale1,16*exp(P.R(3)));
        
    case 3
        % NOISE
        %------------------------------------------------------------------
        rng default;
        mu    = exp(P.R(1));              % mean amplitude
        hfn   = randn(length(pst),1) + (sqrt(-1)*randn(length(pst),1)*1/32);
        drive = hfn*mu;   % amplitude (constant) over time
        drive = mu.*anoise(w,2*length(pst));
        drive = drive(1:length(pst));
        
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

        
        [Pks,PksAmps]=findpeaks(Gu,w);
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
    [y{c},w,s{c},g{c},t{c},layers{c},noise{c},firing{c},QD{c},Spike{c},condel{c}] = ...
        dodxdt(pst,f,v,Q,M,dt,w,drive,Kx,U,method,solvefp);

end

other.noise = noise;
other.firing = firing;
other.QD = QD;
other.Spike = Spike;
other.drive = drive;
other.condel=condel;
end

function [y,w,s,g,t,layers,noise,firing,QD,spike,condel] = ...
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

if isfield(M,'IntMethod')
    IntMethod = M.IntMethod;
end

% Prerequisits for integration with delays
if WithDelays == 2 || WithDelays == 5 || WithDelays == 20 || WithDelays == 21 ...
        || WithDelays == 22 || WithDelays == 8 || WithDelays == 24 || WithDelays == 101 
    
    % Call the function 
    [fx, dfdx,D] = f(M.x,4,P,M);
    
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
    Q       = (spm_expm(dt*D*dfdx/N) - speye(n,n))*spm_inv(dfdx);% system delay operator, Q
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

if npp == 8
    del = exp(P.ID).*[2 1/4 1/2 8 1/2 4 2 2]/2.4;
    del = repmat(del,[1 nk]);
    del=1./del;
    if ns > 1
        if ~isfield(P,'delay')
            del = (spm_vec(repmat(del,[ns 1])))';
        else
        end
    end
    condel=del;
end
if npp == 4 % cmc
    del = exp(P.ID).*[2 1/4 1/2 4]/2.4;
    del = repmat(del,[1 nk]);
    del=1./del;
    if ns > 1
        if ~isfield(P,'delay')
            del = (spm_vec(repmat(del,[ns 1])))';
        else
        end
    end
    condel=del;
end
 
% initial firing rate
Curfire = zeros(size(M.x,2),1)';
firings = [];
DoSpecResp = 1;

% Initialise states [v] and period [t] for integration loop
%--------------------------------------------------------------------------
M.pst = t;
v     = spm_vec(M.x);

if WithDelays == 21
    % pre reqs. for int with bilinear jac-get Jacobian and its derivatives
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

% if solvefp
%    v    = spm_vec(atcm.fun.solvefixedpoint(P,M,drive(1)));
% end

% Frequency steps: dw
dw = 1./(w(2)-w(1));


% % do a burn-in - 
% for i = 1:1200
%     burny(i,:) = spm_vec(v + del'.*(Q*f(spm_unvec(v,M.x),drive(i,:),P,M)));
%     v = spm_vec(burny(i,:));
% end
% v = spm_robust_average(burny,1)';

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
        % integrate: x(t) = T*exp(D*t)*iT*x0        
        for i = 1:length(t)
            y(:,i) =  T*exp(D*(t(i)./1000))*iT*x0;
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
                if drive(i) == 0
                    % baseline
                    v = v;
                else
                    
                    for j = 1:N
                        
                        v = v + del'.*(Q*f(spm_unvec(v,M.x),drive(i,:),P,M));           
                        
                        %[f0,dfdx,D] = f(spm_unvec(v,M.x),drive(i,:),P,M);
                        %Q       = (spm_expm(dt*D*dfdx/N) - speye(n,n))*spm_inv(dfdx);
                        
                        %v = v + del'.*Q*f0;
                        
                        % Ozaki 1992 numerical method:
                        % "dx(t) = (expm(dfdx*t) - I)*inv(dfdx)*f"
                        %[fx,dfdx] = f(v,drive(i),P,M);
                        %v = v + spm_dx(D*dfdx,D*fx,dt);
                    end   
                    
                end
                
                % Expansion point - i.e. deviation around fixed point
                if ~IsStochastic
                    y(:,i) = v - spm_vec(M.x);                    
                else
                    y(:,i) = v - spm_vec(M.x) + rand(size(spm_vec(M.x)));
                end
                
                
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
            fired     = find(squeeze(V(:,:,1)) >= VR);
            firings   = [firings; [i+0*fired',fired'] ];
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


[y,s,g,noise,layers] = spectral_response(P,M,y,w,npp,nk,ns,t,nf,timeseries,dt,dfdx);
t = drive;
end

function [y,s,g,noise,layers]=spectral_response(P,M,y,w,npp,nk,ns,t,nf,timeseries,dt,dfdx)

% Spectral Response Options
%--------------------------------------------------------------------------
%DoSmoothing    = 1; % smooth the fourier series with a 10% Loess window
DoHilbert      = 0; % take the absolute (magnitude) of the hilbert envelope
Bandpassfilter = 0; % band pass filter to [w(1)-1) : w]
DoDCT          = 0; % discrete cosine transform series before fft
IncDCS         = 1; % include semi-stochastic neuronal fluctuations       x       % ON FOR DEXPRO
DoHamming      = 1; %(1)% Cosine Hann window the model spectrum      
HamNoise       = 0; % Hann noise components - if data was BPF, exponential 
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

% Neuronal innovations: a multiplier on the model signal
%--------------------------------------------------------------------------
for i = 1:size(P.a,2)
    Gu(:,i) = exp(P.a(1,i))*(w.^0);                  % P.a = constant
end

% Spectrum of channel noise (non-specific): added to spectrum
%--------------------------------------------------------------------------
Gn = exp(P.b(1,i) )+w.^(-exp(P.b(2,1))); 

% Spectrum of channel noise (specific): added to spectrum
%--------------------------------------------------------------------------
for i = 1:size(P.c,2)
    Gs(:,i) = exp(P.c(1,i) )+w.^(-exp(P.c(2,1)));    % P.c = expone
end
if HamNoise
    % Hamming to taper the edges 
    % H  = (1 - cos(2*pi*[1:nf]'/(nf + 1)))/2;
    % optimisation struggles with edge effects
    %----------------------------------------------------------------------
    warning off ;     % dont warn of integer operands
    H  = kaiser(nf,2.5);
    Gs = (Gs.*H);
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
burn = findthenearest(300,t); 

if isfield(M,'burnin')
    burn = findthenearest(M.burnin,t); 
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

% Select the frequency transfer function - or default smoothed fft(x)
%--------------------------------------------------------------------------
fmethod = 'none'; % none, svd, dmd, timefreq ...
if isfield(M,'fmethod')
    fmethod = M.fmethod;
end

if M.timefreq
    % force, non-negotiable...
    fmethod = 'timefreq';
end

% Whe using DMD remove the ~principal component
try
    switch fmethod            % UNCOMMENT FOR KET PAPER
        case {'dmd' };
            for ins = 1:ns
                % remove principal [dominant] eigenmode(s) from states series
                x = squeeze( y(ins,:,:,:) );
                x = reshape( x , [npp*nk, length(t)] );
                [u0,s0,v0] = spm_svd(((x)));
                %ipc = find( (cumsum(diag(full(s0)))./sum(diag(full(s0))) > .8) );
                p1 = u0(:,1)*s0(1,1)*v0(:,1)';
                p2 = u0(:,2)*s0(2,2)*v0(:,2)';
                x = x - (p1 + p2);
                y(ins,:,:,:) = reshape(x,[npp,nk,length(t)]);
            end
            timeseries=y;
    end
end

dw = 1./(w(2)-w(1));

% This is the main loop over regions to calculate the spectrum
%--------------------------------------------------------------------------
for ins = 1:ns
        
    % extract time series of all states from this region
    %----------------------------------------------------------------------
    yx = reshape( squeeze(y(ins,:,:,:)), [npp*nk,length(t)]); 

    % parameter controlled suppression of cells: P.Jpop(i)
    if isfield(P,'Jpop') && length(P.Jpop) == 8
        for pij = 1:8
            yx(pij,:) = yx(pij,:)*P.Jpop(pij);
        end
    end
        
    % use dynamic mode decomposition (DMD) to reduce the oscillatory
    % complexity to a simplistic model: a set of discrete spatio-temporal modes
    %----------------------------------------------------------------------
    warning off;              % avoid rank deficiency warning
    switch lower(fmethod)
        case {'kpca' 'dmd' 'instantaneous' 'svd' 'none' 'glm' 'fooof' 'timefreq' 'eig'}
            switch fmethod
                
                case 'timefreq'
                    
                    if ~isfield(M,'nw')
                        M.nw = 30;
                    end
                    
                    tfmat = 0;
                    for ij = 1:length(Ji)
                        % load the virtual sensor and run a timefrequency analysis
                        cfg.baseline = 'relchange';
                        cfg.sampletimes = double(M.pst);
                        cfg.fsample = 1./dt;
                        cfg.filterorder = 4;
                        FoI = linspace(w(1),w(end),M.nw);

                        %MatDat = real(J(:)'*yx);
                        MatDat = real(yx(Ji(ij),:));
                        
                        tf{i} = atcm.fun.bert_singlechannel([MatDat],cfg,FoI,[-1 0]);
                        y = double(tf{i}.agram);
                        y = double(atcm.fun.HighResMeanFilt(y,1,2));
                        tfmat = tfmat + ( Ji(ij)*y );
                    end
                    
                    y = exp(P.L)*tfmat;
                    s = yx;
                    g = ts;
                    noise = [];
                    layers.iweighted = yx;
                    return;                
                
                case 'kpca'
                    
                    x = yx(1:8,:);
                    kpca = KernelPca(x', 'gaussian', 'gamma', 2.5, 'AutoScale', true);
                    M0 = length(find(Ji));
                    Eigenvectors = project(kpca, x', M0)';
                    
                    
                case 'dmd'
                [Eigenvalues, Eigenvectors, ModeAmplitudes, ModeFrequencies, ...
                    GrowthRates, POD_Mode_Energies] = atcm.fun.dmd((yx), length(Ji), dt);
%                 % Manual DMD:
%                 X1 = yx(:,1:end-1);
%                 X2 = yx(:,2:end);
%                 [U,S,V] = spm_svd(X1);
%                 r = length(Ji);
%                 U = U(:,1:r);
%                 S = S(1:r,1:r);
%                 V = V(:,1:r);
%                 Atilde = full(U')*X2*full(V)*pinv(full(S));
%                 [W,eigs] = eig(Atilde);
%                 Phi = X2*V*pinv(full(S))*W;
%                 Eigenvectors = Phi'*yx;
                case 'instantaneous'                    
                    yx(1,:) = yx(1,:)*0.2;
                    yx(2,:) = yx(2,:)*0.8;
                    yx(3,:) = yx(3,:)*0.1;
                    yx(4,:) = yx(4,:)*0.2;
                    yx(5,:) = yx(5,:)*0.1;
                    yx(6,:) = yx(6,:)*0.2;
                    yx(7,:) = yx(7,:)*0.05;
                    yx(8,:) = yx(8,:)*0.1;
                    [Eigenvectors,s,v] = svd(yx);
                    Eigenvectors = Eigenvectors'*yx;       
                case 'svd'
                   % [Eigenvectors,s,v] = svd(cov(real(yx)'));
                    [Eigenvectors,s,v] = svd(real(yx));
                    Eigenvectors = Eigenvectors'*yx; 
                case {'none' 'fooof' 'timefreq' 'eig'}
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
                                         
                        
                        case {'timefreq'}
                            data = Eigenvectors(Ji(ij),burn:end);
                        
                        case 'eig'
                                J = dfdx;
                                [v,s] = eig(full(J),'nobalance');
                                s     = diag(s);
                                s     = 1j*imag(s) + real(s) - exp(real(s));
                                
                                [vv,uu,ss]=svd(Eigenvectors*Eigenvectors');
                                dvdu  = pinv(vv);

                                for k = 1:length(s)
                                    Sk       = 1./(1j*2*pi*w - s(k));
                                    %S(k,:) = dvdu(k,:)*Sk;
                                    SS(k,:) =Sk;
                                end
                                spec = dvdu*SS;
                                Pf = exp(P.J(Ji(ij))) * spec(Ji(ij),:)';
                            
                        case {'none','dmd','svd' 'kpca'}
                            % just a smoothed fft of the (contributing)states
                            switch fmethod
                                case 'none';this = Eigenvectors(Ji(ij),burn:end);
                                case 'dmd' ;this = y0(burn:end);
                                case 'svd' ;this = Eigenvectors(ij,burn:end);
                                case 'kpca';this = y0(burn:end);
                            end
                            
%                             switch fmethod
%                                 case 'none'
%                                     if ij == 1
%                                         this = ts(:);
%                                     else
%                                         continue;
%                                     end
%                             end
                            
                            ncompe = 100;
                            if isfield(M,'ncompe')
                                ncompe = M.ncompe;
                            end
                            
                            % An AR routine
%                            this2 = spm_vec(this)*spm_vec(this)';
%                            R = fliplr(eye(length(this2)));
%                            AC = (1/2)*(this2 + R*this2'*R);
%                            [V,s] = eig(AC);
%                            fq = maxpoints((real(log(imag(diag(s))))) , 10 );
%                            for k = 1:length(fq)
%                                Sk(k,:) = 1./(1j*2*pi*w - fq(k));
%                            end
%                            Pf = mean(Sk,1)';

                            % pad & filter
                       %     [thispad,It] = atcm.fun.padtimeseries(this);
                       %     thispad = atcm.fun.bandpassfilter(thispad',1./dt,[w(1)-1 w(end)]);
                       %     this = thispad(It);
%                             
%                             % original
%                             PfOrig  = atcm.fun.Afft(this,dw/dt,w);
%                             %Pf = pyulear(this,12,w,dw./dt);%.*Hz.^2;
%                             nwg = 4;%*exp(P.psmooth);
%                             w0 = 1 + (nwg*( w./w(end)));
%                             PfOrig = PfOrig(:).*w0(:);%.*Hz(:);
%                             
%                             [padp0,indp0] = atcm.fun.padtimeseries(PfOrig);
%                             Pfs0 = atcm.fun.tsmovavg(padp0','t',8);
%                             Pfs0 = Pfs0(indp0);
%                             PfOrig=Pfs0(:);
%                             
%                             NC = 30;
%                             
%                             clear Pf
%                             
%                             % decompose timeseries using SSA
%                             sthis = atcm.fun.assa(this,NC);
%                             
%                             for ikw = 1:NC
%                                 [Pf,Hz]  = atcm.fun.Afft(sthis(:,ikw)',dw/dt,w);
%                                 nwg = 4;%*exp(P.psmooth);
%                                 w0 = 1 + (nwg*( w./w(end)));
%                                 PfPf(ikw,:) = Pf(:).*w0(:);%.*Hz(:);
%                             end
%                             
%                             DoEnv = 1;
%                             if isfield(M,'DoEnv')
%                                 DoEnv=M.DoEnv;
%                             end
%                             
%                             if DoEnv
%                                 for ikw = 1:NC
%                                     [padp,indp] = atcm.fun.padtimeseries(PfPf(ikw,:));
%                                     Pfs = atcm.fun.tsmovavg(padp','t',8);
%                                     Pf1 = Pfs(indp);
%                                     Pf(ikw,:)=Pf1(:);
%                                 end
%                             else
%                                 Pf = PfPf;
%                             end
%                             
%                             
%                             pcx = corr(Pf',PfOrig).^2;
%                             %[~,I]=atcm.fun.maxpoints(pcx,4);
%                             [~,I]=sort(pcx,'descend');
%                             these = atcm.fun.findthenearest(cumsum(pcx(I))./sum(pcx),.9);
%                             I = I(1:these);
%                             
%                             Pf = sum(Pf([I],:),1)';
                            
                            
                            
%                             
%                             % only explain away error
%                             %if ij == 1
%                                 yy = spm_vec(M.y);
%                             %else
%                             %    yy = spm_vec(yy) - Pfx(:);
%                             %end
%                             
%                             %AllPfComp{ij} = Pf;
%                             
%                             warning off;
%                             Bb = atcm.fun.lsqnonneg(Pf',yy);
%                             warning on;
%                             
%                             Pf = spm_vec(Bb'*Pf);
                            
                            %try   Pfx = Pf + Pfx;
                            %catch;Pfx = Pf;
                            %end
   

                            if ncompe > 0
                               [Pf,Hz,Pfmean]  = atcm.fun.AfftSmooth(this,dw/dt,w,ncompe);                             
                               Pfmean = squeeze(Pfmean);
                               Pf = spm_vec(max(Pfmean'));
                               
                               nwg = 4*exp(P.psmooth);
                                w0 = 1 + (nwg*( w./w(end)));
                                Pf = Pf(:).*w0(:);
                            else
                                [Pf,Hz]  = atcm.fun.Afft(this,dw/dt,w);
                                %Pf = pyulear(this,12,w,dw./dt);%.*Hz.^2;
                                nwg = 4;%*exp(P.psmooth);
                                w0 = 1 + (nwg*( w./w(end)));
                               % w0 = linspace(1,2,length(w));
                                Pf = Pf(:).*Hz(:);%w0(:);%.*Hz(:);
                            end
                                                       
                            DoEnv = 1;
                            if isfield(M,'DoEnv')
                                DoEnv=M.DoEnv;
                            end
                            
                            if DoEnv
                                
                               [padp,indp] = atcm.fun.padtimeseries(Pf);
                               Pfs = atcm.fun.tsmovavg(padp','t',8);
                               Pf = Pfs(indp);
                               Pf=Pf(:);
                               
                            end
                            
                            % store 
                            Pf0(ins,ij,:) = Pf;
                                                        
                        case {'fooof'}
                            
                            fq = w(1):(1/8):w(end);
                            [FT,~]  = atcm.fun.Afft(y0(burn:end),1/dt,fq);
                            res = fooof(fq, FT, w, struct(), 1);
                            
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

                J  = full(exp(P.J));

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
    case {'svd','dmd'};
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

% Incorporate noise components for auto (Gs) and cross (Gn) spectra
%----------------------------------------------------------------------        
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

%clear Pf;

% If M.y contains the empirical data, fit it as a GLM of the contirbuting
% populations
if isfield(M,'y')
    for ins = 1:ns
        dat = squeeze(layers.iweighted(ins,:,:))';
        yy  = squeeze(M.y{1}(:,ins,ins));
        if all( (size(dat,1)==1) && size(dat,2)==length(w) )
            dat = dat';
        end
        linmod = 1;
        if isfield(M,'linmod')
            linmod = M.linmod;
        end
        addnoise = 1;
        if linmod == 1
            if isfield(M,'envonly') && M.envonly == 1
                Mm = dat';
                b = exp(P.J(Ji));
                %b  = pinv(Mm*Mm')*Mm*yy;
                if any(size(Mm)==length(b))
                    Pf(:,ins,ins) = b(:)'*Mm; 
                    Pf(:,ins,ins) = Pf(:,ins,ins) * exp(P.L(ins));
                else
                    Pf(:,ins,ins) = Mm;
                end
                
            elseif isfield(M,'envonly') && M.envonly == 2
                addnoise = 0;
                Mm = dat';
                Mm = [Mm; Gu'; Gn]; % put the noise in the glm, not added aftwrwards
                b  = pinv(Mm*Mm')*Mm*yy;
                Pf(:,ins,ins) = b'*Mm; 
                Pf(:,ins,ins) = Pf(:,ins,ins) * exp(P.L(ins));
            end
        end
        
        if isfield(M,'EnvLFP') && M.EnvLFP
            if isfield(M,'LFPsmooth') && ~isempty(M.LFPsmooth)
                smthk = round(M.LFPsmooth);
            else
                smthk = 2;
            end
                        
            if isfield(M,'tsmovavg')
                do_movavg = M.tsmovavg;
            else
                do_movavg = 0;
            end
            
            if ~do_movavg
                i0 = Pf(:,ins,ins);
                yy ;
                usesmoothkernels = 0;
                if isfield(M,'usesmoothkernels') && M.usesmoothkernels
                    usesmoothkernels = 1;
                end
                
                if ~usesmoothkernels
                                        
%                     % Compute singular spectrum analysis [ssa], fir comps
%                     %------------------------------------------------------
                    X = Pf(:,ins,ins);
                    
                    %[Xp,Xpi] = atcm.fun.padtimeseries(X);
                   % Xp = full(atcm.fun.HighResMeanFilt(Xp,1,smthk));
                   % X = Xp(Xpi);
                   % Pf(:,ins,ins)=X;
                    
                   RC = atcm.fun.assa(X,30);
                   
                    pc = RC;
                    %Bw = pinv(pc'*pc)*pc'*yy;
                                        
                    %pcx = corr(pc,Pf(:,ins,ins)).^2;
                    weight = M.FS(M.y{:});
                    %weight = weight./sum(weight);
                    %pcx = atcm.fun.wcor([pc Pf(:,ins,ins)],w'./sum(w)).^2;
                    pcx = atcm.fun.wcor([pc Pf(:,ins,ins)],weight).^2;
                    pcx = pcx(1:end-1,end);
                    %[~,I]=atcm.fun.maxpoints(pcx,4);
                    [~,I]=sort(pcx,'descend');
                    these = atcm.fun.findthenearest(cumsum(pcx(I))./sum(pcx),.4);
                    I = I(1:these);
                    %fprintf('%d/%d\n',these,length(pcx));
                    
                    Pf(:,ins,ins) = exp(P.L(ins))* sum(pc(:,[I]),2);
                    
                    %[padpf,ipf] = atcm.fun.padtimeseries(Pf(:,ins,ins));
                    %padpf = full(atcm.fun.HighResMeanFilt(padpf,1,smthk));
                    
                    %Pf(:,ins,ins) = exp(P.L(ins))*padpf(ipf);
                    
                    %warning off;
                    %Bw = lsqnonneg(pc,yy); % pos constr LSQGLM
                    %warning on;
                                        
                    %Pf(:,ins,ins) = exp(P.L(ins))*(pc*Bw);
                    
                   % warning off;
                   % Bw = lsqnonneg(cat(1,AllPfComp{:})',yy); % pos constr LSQGLM
                   % warning on;
                                        
                   % Pf(:,ins,ins) = exp(P.L(ins))*(cat(1,AllPfComp{:})'*Bw);

                 %   Pf(:,ins,ins) = exp(P.L(ins)) * Pf(:,ins,ins);
                    
                    %c = RC;
                    %X = {pc Bw};
                    c=[];
                    X=[];
                    meanpower={X c};
                    
%                     m = fit(w.',Pf(:,ins,ins),'fourier8');
%                     x = Pf(:,ins,ins);
%                     p = m.w; %Dirichlet's condition 
%                     
%                     % rebuild the cossines constituting the fourier series
%                     for j = 1:8
%                         c(j,:) = m.(['a' num2str(j)])*cos(j*w*p)+m.(['b' num2str(j)])*sin(j*w*p);
%                     end
%                     warning off;
%                     X = atcm.fun.lsqnonneg(c',yy); % pos constr LSQGLM
%                     %X = pinv(c*c')*c*yy;
%                     %Pf(:,ins,ins) = exp(P.L(ins)) * m(w);
%                     warning on;
%                     Pf(:,ins,ins) = exp(P.L(ins)) * (X'*c);
                    %Pf = mod(w);
                    
                    
                    %X=[];
                    %c=[];
                    
%                     
% %                     % optimise smoothing function by picking best from an
% %                     % iteratively smoothed version
%                      [padpf,ipf] = atcm.fun.padtimeseries(Pf(:,ins,ins));
%                      smoothpf = full(atcm.fun.HighResMeanFilt(padpf,1,smthk));%12
% %                     smoothpf = atcm.fun.tsmovavg(smoothpf','e',smthk);
% %                     smoothpf = abs(hilbert(smoothpf));
%                     i1 = smoothpf(ipf);
%                     
%                     for ik = 1:length(i0)
%                         opts(ik,:) = linspace(i0(ik),i1(ik),4);
%                         
%                         %[pado,sec] = atcm.fun.padtimeseries(opts(ik,:));
%                         %pado = abs(hilbert(pado));
%                         %opts(ik,:) = pado(sec);
%                         
%                     end
% % %                     
% % %                     
% % %                     %b = pinv(opts'*opts)*opts'*yy;
% % %                     %out(:,ins,ins) = (opts*b)+Gn(:);
% %                     
% %                     %X = lsqnonneg(opts,yy); % pos constr LSQGLM
% %                     %out(:,ins,ins) = exp(P.L(ins)) * (X'*opts');
% %                    % Pf(:,ins,ins)=out(:,ins,ins);
%                     for ik = 1:length(i0)
%                         this = yy(ik);
%                         opi  = opts(ik,:);
%                         ind  = atcm.fun.findthenearest(this,opi);
%                         ind  = ind(1);
%                         
%                         out(ik,ins,ins)=opi(ind);
%                     end
%                      Pf(:,ins,ins) = exp(P.L(ins))*out(:,ins,ins);
                   %  [padvec,yi] = atcm.fun.padtimeseries(Pf(:,ins,ins));
                   %  tsmth = full(atcm.fun.HighResMeanFilt(padvec,1,smthk));
                %     tsmth = atcm.fun.tsmovavg(tsmth','e',smthk);
                %     %tsmth = full(atcm.fun.tsmovavg(padvec','t',smthk));
                    % Pf(:,ins,ins) = exp(P.L(ins))*tsmth(yi);
                    
                    %Pf(:,ins,ins) = full(atcm.fun.HighResMeanFilt(Pf(:,ins,ins),1,smthk));
                  %  Pf(:,ins,ins) = exp(P.L(ins))*atcm.fun.aenvelope(Pf(:,ins,ins),20);
                else
                    [APf,GL] = AGenQ(Pf); % AGenQ gens smoothing vectors & GraphLap       
                    APf = APf;%.*GL;
                    i0  = Pf(:,ins,ins);
                    for ik = 1:length(i0)
                        this = yy(ik);
                        opi  = APf(ik,:).*Pf;
                        ind  = atcm.fun.findthenearest(this,opi);
                        ind  = ind(1);
                        out(ik,ins,ins)=opi(ind);
                    end
                    Pf(:,ins,ins) = out(:,ins,ins);
                    Pf(:,ins,ins) = full(atcm.fun.HighResMeanFilt(Pf(:,ins,ins),1,smthk));
                    Pf(:,ins,ins) = exp(P.L(ins))*atcm.fun.aenvelope(Pf(:,ins,ins),35);
                end
            else
                Ppf{1} = atcm.fun.tsmovavg(Pf(:,ins,ins)','e',2)';
                Ppf{2} = atcm.fun.tsmovavg(Pf(:,ins,ins)','e',4)';
                Ppf{3} = atcm.fun.tsmovavg(Pf(:,ins,ins)','e',8)';
                Ppf{4} = atcm.fun.tsmovavg(Pf(:,ins,ins)','e',12)';
                
                opt = [Pf(:,ins,ins) cat(2,Ppf{:})];
                for iy = 1:length(yy)
                    [~,I]=min(( yy(iy) - opt(iy,:)).^2);
                    I = I(1);
                    Pf(iy,ins,ins) = opt(iy,I);
                end
            end
        end
        %addnoise=0;
        if addnoise
            % Multiply in the semi-stochastic neuronal fluctuations
            for i = 1:length(Hz)
                Pf(i,:,:) = sq(Pf(i,:,:))*diag(Gu(i,ins));
            end
        end
        %Pf(:,ins,ins) = smooth( squeeze(Pf(:,ins,ins)) , 4*exp(P.psmooth(1)) ,'moving' );
    end

    for inx = 1:ns
        for iny = 1:ns
            if inx ~= iny
                Pf(:,inx,iny) = squeeze( Pf(:,inx,inx) ) .* ...
                                   conj( Pf(:,iny,iny) ) ;
            end
        end
    end
    
    % RE-Incorporate noise components for auto (Gs) and cross (Gn) spectra
    %----------------------------------------------------------------------
    %if ns > 1
    addnoise=1;
    if addnoise
        for i = 1:ns
            for j = 1:ns
                % Autospectral noise / innovations
                Pf(:,i,j) = Pf(:,i,j) + Gn(:);
                if j ~= i
                    % Cross spectral noise / innovations
                   % Pf(:,j,i) = Pf(:,j,i) .* Gs(:,i);
                   % Pf(:,i,j) = Pf(:,j,i);
                end
            end
        end
     end
    
end
%DoHamming=0;
if DoHamming
    for i = 1:ns
        for j = 1:ns
            H  = (1 - cos(2*pi*[1:nf]'/(nf + 1)))/2;
            H  = kaiser(nf,2.5);
            Pf(:,i,j) = Pf(:,i,j).*H;
        end
    end
end

% returns for this trial - {g}
%--------------------------------------------------------------------------
y = Pf;
s = timeseries;
g = ts;
try
    noise.aperiodic = meanpower;
end

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
