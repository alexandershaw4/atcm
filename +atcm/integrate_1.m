function [y,w,s,g,t,pst,layers,other] = integrate_1(P,M,U,varargin)
% Numerical integration and spectral response of a neural mass model.
% This version, which can handle networks of connected models,
% and with different trial types. This version (atcm.integrate_1) has a number 
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
% and spectral response from a neural mass, using numerical methods.
%
% The optimisation problem becomes highly nonlinear so i suggest using AO.m
% optimisation (https://github.com/alexandershaw4/aoptim).
%
% This routine is more basic than the deault DCM transfer functions because, 
% rather than taking an fft of the system's kernels (eigenvectors of the 
% Jacobian), it uses the  approach of integrating ( dx/dt ) using a Euler 
% or Runge-Kutta scheme over ~3s to first generate a states time series 
% (voltages, currents). Finally, weighting & an fft of  the resultant 
% timeseries provide the spectral response.
% This is nice because one can access both the continuous and discrete 
% time responses of the system. Channel noise (exponetial decay functions)
% are computed and returned.
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
%--------------------------------------------------------------------------
% All options are selected by specifying fields in the input 'M' structure.
%
% Selecting Input Type:
% -------------------------------------------------------------------------
% By default, the input is a constant (D.C).
% Set M.InputType = 0 ... constant
%                 = 1 ... a sine wave oscillation
%                 = 2 ... an ERP (Guassian) bump
%                 = 3 ... high frequency noise
%
% Select Sample Period:
% -------------------------------------------------------------------------
% Default = 3s @ 600Hz.
% Set M.sim.pst = vector of sample times
%     M.sim.dt  = update step (1./fs)
%
% Selecting Numerical Integration Method:
% -------------------------------------------------------------------------
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
%                 = 44 ... Newton-Cotes with delays
%                 = 21 ... integration with bilinear jacobian
%                 = 23 ... full jacobian integration
%                 = 24 ... stochastic equation integration
%                 = 8  ... 8th-order Runge-Kutta w/ delays
%                 = 45 ... 4th order Runge-Kutta w/ delays
%
% Time series decomposition 
% -------------------------------------------------------------------------
% Set M.decompose = 'ssa'   - singular spectrum analysis algorithm
%                 = 'fourier' - fit a fourier series
%
%
% Selecting Spectral (transfer) Method
% -------------------------------------------------------------------------
% Set M.fmethod = 'dmd' ... use dynamic mode decomposition to identify 
%                           frequency modes in the intgrated timeseries
%               = 'none' ... just use an fft (or smooth fft) of the data
%
% Note - when using 'none', optionally specify how many smoothing windows
% to use in M.smth, default = 30
%
% Spectral smoothing
% -------------------------------------------------------------------------
% The fft(x) for each state is computed using a sliding window-smoothed
% fourier transfer (atcm.fun.AfftSmooth) on a DCT filtered (retaining 99%)
% version of the timeseries. The resulting spectrum is fit with a contrained
% exponential (atcm.fun.c_oof) or some combinatorial multivariate distribution
% model (using Gaussian, Cauchy, Laplace and Gamma dists) 
% 
% Other options:
% -------------------------------------------------------------------------
% M.IncDCS = flag to include a discrete cosine set (a semi-stochastic set
% of neuronal fluctuations in frequency space).
%
% Also required: SPM12 w/ DCM,  
%
% Dr Alexander Shaw | 2020 | alexandershaw4[@]gmail.com


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
    dt    = dt*exp(P.dt);
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

% Select input type: 0 = constant (DC), 1 = oscilltion (sine), 2 = ERP bump
% 3 = noise, 4 = 2 oscillations, 5 = ifft of parameterised basis set
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
        mu    = .2*exp(P.R(1));                      % mean amplitude
        mf    = 100*exp(P.R(2));                      % frequency
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
        %rng default;
        mu    = exp(P.R(1));              % mean amplitude
        hfn   = randn(length(pst),1) + (sqrt(-1)*randn(length(pst),1)*1/32);
        drive = hfn*mu;   % amplitude (constant) over time
        drive = .25*drive(1:length(pst));
        
    case 4
        % TWO oscillatory inputs...
        %------------------------------------------------------------------
        mu1   = 0.5*exp(P.R(1));                % mean amplitude
        mu2   = 1.0*exp(P.R(2));
        mf1   = 100*exp(P.R(3));                  % frequency
        mf2   = 1*exp(P.R(4));
        
        drive(:,2) = (mu1 * sin(2*pi*mf1*(pst/1000)) ) ;  
        drive(:,1) = (mu2 * sin(2*pi*mf2*(pst/1000)) );
        
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
        
    case 6

%         % ifft a parameterised DCT set: neuronal fluctuations
%         for i = 1:size(P.a,2)
%             Gu(:,i) = exp(P.a(1,i))*(w.^0);        % P.a = constant
%         end
%         nf = length(w);
%         nd = size(P.d,1);
%         X  = spm_dctmtx(nf,nd + 1);
%         Mu = exp(X(:,2:end)*P.d);
%         if size(Mu,2) == 1, Mu = Mu*ones(1,1); end
%         Gu = Gu.*Mu;
%         Gu = exp(P.a(1))*Mu;
% 
%         [Pks,PksAmps]=findpeaks(real(Gu),w);
%         for i = 1:length(Pks)
%             gx(i,:) = PksAmps(i) * sin(2*pi*Pks(i)*pst);
%         end
%         drive(:,2) = sum(gx,1)';
        
        % For oscillatory inputs...
        %------------------------------------------------------------------
        mu    = .2*exp(P.R(1));                      % mean amplitude
        mf    = 10*exp(P.R(2));                      % frequency
        drive(:,2) = mu * ( (sin(2*pi*mf*(pst/1000))) );%...
        %   + sin(2*pi*(10*exp(P.R(3)))*(pst/1000)) );
                  

        % and the erp
                
        % For ERP inputs...
        %------------------------------------------------------------------
        delay  = 60 * exp(P.R(1));             % bump
        scale1 = 8  * exp(P.R(2));
        drive(:,1)  = atcm.fun.makef(pst,delay,scale1,16*exp(P.R(3)));
        drive(1,1)=0;
        
        sust = (max(drive(:,1)))*.75;
        intcpt = atcm.fun.findthenearest(drive(:,1),sust);
        drive(intcpt:end,1) = sust;%*wave(intcpt:end);

end


if ~isfield(M,'timefreq')
    M.timefreq = 0;
end

% expansion (fixed) point: trial & parameter effects are deviations from here
%--------------------------------------------------------------------------
f    = spm_funcheck(M.f); 

% solve for a fixed point, or not
if solvefp; 
    %x    = atcm.fun.solvefixedpoint(P,M);
    x = atcm.fun.solvefixedpoint(P,M,[],-70);
else ;      x    = x;
end

M.x  = x;
v    = spm_vec(x);
NoFX = 0;

% flag no modulatory effects in this model
if isempty(U)
    U.X  = 1;
    NoFX = 1; 
end

% Integration and spectral response for this trial (c)
%--------------------------------------------------------------------------
for  c = 1:size(U.X,1)
    
    % generate condition-specific parameter structure
    %----------------------------------------------------------------------
    if ~NoFX; Q  = spm_gen_Q(P,U.X(c,:));
    else      Q  = P;
    end
        
    % integration, spectral response, firing etc. (subfunctions)
    %----------------------------------------------------------------------
    [y{c},w,s{c},g{c},t{c},layers{c},noise{c},firing{c},QD{c},Spike{c},condel{c},series{c}] = ...
        dodxdt(pst,f,v,Q,M,dt,w,drive,Kx,U,method,solvefp,c);

end

% outputs
%--------------------------------------------------------------------------
other.noise = noise;
other.firing = firing;
other.QD = QD;
other.Spike = Spike;
other.drive = drive;
other.condel=condel;
other.series = series;
other.dt = dt;
other.Fs = Fs;
end

function [y,w,s,g,t,layers,noise,firing,QD,spike,condel,series] = ...
                            dodxdt(t,f,v,P,M,dt,w,drive,Kx,U,method,solvefp,ci)
% Numerical integration, signal weighting and FFT of timeseries with
% spline interpolation and smoothing


% Choose either the DCM kernels scheme, or another dx = f(x) update scheme
%--------------------------------------------------------------------------
if isempty(method)
    method = 2;
end

warning off;
IntMethod   =  ' ';
WithDelays  = method;  
[ns,npp,nk] = size(M.x);

if isfield(M,'IntMethod')
    IntMethod = M.IntMethod;
end

% Prerequisits for integration with delay operators
%--------------------------------------------------------------------------
[fx, dfdx,D] = f(M.x,[4 4],P,M);
%dFdx = dfdx;

  OPT.tol = 1e-6*norm((dfdx),'inf');
  if OPT.tol == 0
       OPT.tol = 1;
  end
 p     = abs(eigs(dfdx,1,'SR',OPT));
 N     = ceil(max(1,dt*p*2));

% if isfield(P,'ID')
%     %del = exp(P.ID).*[1 1 1 1 1 1 1 80];
%     
%     id = exp(P.ID).*[1 1 1 1 1 1 1 1];
%     %id = exp(P.ID).*[1 1/2 1/2 1 1 1 1 1]*2;
%     
%     %del = exp(P.ID).*[1 1 1 1 1 1 1 80];
%     
%     CT = 8*exp(P.CT); %60;
%     TC = 3*exp(P.TC); %20;
%     
%     tc = [1 1 1 1 1 1 TC TC];
%     ct = [CT CT CT CT CT CT 1 1];
%     
%     del = zeros(8,8);
%     
%     del(1:6,[7 8]) = TC;
%     del([7 8],1:6) = CT;
%    
%     del = del + diag(id);
%     
%     %del = -del / 1000;
%     
%     
% else
%     del = ones(1,8);
% end
% 
% %del = repmat(del,[1 nk]);
% del = repmat(del,[nk nk]);

%del = 1./del;

% Convert delay matrix to update-rate matrix - for the ensuing update algo:
%--------------------------------------------------------------------------
%
%   x(t+1) = x(t) + dt * [ Q*dfdx-I * f(x,...) ]
% 
% where rate matrix Q[k] == 1 if state k has zero delay, is just less than
% 1 if small delay and is a very small number if Q[k] is a large delay.

% Q  = diag( 1./(1 - (dt*del))-dt ); D = Q; % <-- NOTE!
% U  = spm_svd(dfdx,0);
% Q  = U'*(Q*dfdx)*U;

%dfdx = U'*dfdx*U;

%D = ( 1./(1 - (dt*del))-dt ); D = D + dt;
%D = del;
%D=[];
% regularise Jacobian
%dfdx = dfdx + ~eye(length(dfdx))*min(dfdx(dfdx~=0));

%dfdx = full(real(dfdx)) + ~eye(56)*1e-2;
%dfdx(abs(dfdx)<1)=1;

%dfdx = rescale(dfdx,-2,2);

%dfdx = (dfdx + dfdx') /2;

%Q     = (spm_expm(dt*D.*dfdx/N) - speye(56,56))/dfdx;
Q     = (spm_expm(dt*D.*dfdx/N) - speye(56,56))/dfdx;

%Q     = (spm_expm(dt*(D.*dfdx))  )./dfdx;
%Q = (D.*dfdx)*dt;

Q = D*dt;

%U  = spm_svd(dfdx,0);
%Q  = U'*Q*U;

condel = 1;
QD     = Q; 

% firing rate & count when fired (membrane potential passes threshold)
%-------------------------------------------------------------------------
Curfire    = zeros(size(M.x,2),1)';
firings    = [];
DoSpecResp = 1;

% limit delays to make computation of derivatives tractable
try; N = min(N,14); end

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

% % Regularise the Jacobian and delay operators
% if isfield(M,'regularise') && M.regularise
%     lbdmin = min(eigs(dfdx));
%     dfdx = dfdx + 0.5*max(-lbdmin,0)*ones(size(dfdx));
%     lbdmin = min(eigs(cov(Q)));
%     Q = Q + 0.5*max(-lbdmin,0);
% end

switch IntMethod

    case 'kernels'
        
        % Just use the kernels approach?
        %------------------------------------------------------------------
        [K0,K1,K2,H1] = spm_kernels(M,P,length(t),dt);
        y = H1(:,1:end-1)';
        S = [];
        
    case 'ode45'
        % matlab build in ode solver
        ode = @(t,v,P,M,f) spm_vec( f(spm_unvec(v,M.x),drive(1),P,M) );
        
        %ode = @(t,v,P,M,f) spm_vec( Q*f(spm_unvec(v,M.x),drive(1),P,M) );
        
        %[~,y]   = ode113(ode,t/1000,spm_vec(v),drive,P,M,f);
        
        [~,y]   = ode113(ode,[t/1000],spm_vec(v),drive,P,M,f);
        y = y';
        S = [];
        
    case 'linearise'
        % Here we aim to (bi)linearise f, such that:
        % dx = f(x) --> dx = Ax
        % or more precisely,
        % dx = f(x,u) --> dx = Ax + Bu  where A describes the system flow
        %                               and B the input(u)-related dynamics
        %
        % then we use matrices A&B to predict the timeseries of the system
        
        % note: M.x should be at (or close to) a fixed point
        
        % for more info, see Steve Brunton's excellent youtube tutorial:
        % https://www.youtube.com/watch?v=nyqJJdhReiA
        
        % notes:
        % http://faculty.washington.edu/sbrunton/me564/pdf/L06.pdf
        
        % initial point
        x0      = spm_vec(M.x);
        
        % find a fixed point
        if solvefp; xbar    = spm_vec(atcm.fun.solvefixedpoint(P,M));
        else;       xbar    = spm_vec(M.x);
        end
        M.x     = spm_unvec(xbar,M.x);
        
        % compute (states) Jacobian, evaluated at xbar
        dfdx    = spm_diff(M.f,M.x,M.u,P,M,1);
                            
        % compute df/du (inputs)
        for i = 1:length(drive)
            dfdu(:,i)    = spm_diff(M.f,M.x,drive(i),P,M,2);  
        end
        dfdut = dfdu;
        dfdu  = sqrt(dfdu*dfdu');
        
        % eigenvectors and values - dfdx
        [T,D]   = eig(full(real((Q.*dfdx))));
        iT      = pinv(T);
        d       = diag(D);
        
        % eigenvectors and values - dfdu
        [Tu,Du]   = eig(full(real((dfdu))));
        iTu      = pinv(Tu);
        du       = diag(Du);
                
        %in_proj = dfdut*diag(drive);
        
        A = dfdx;
        B = dfdu;
        I = eye(length(dfdx));
        C = exp(P.J);
        
        H = C*inv(I - A)*B; % Numerical Laplace
        
        for Sk = 1:size(H,1)
            Pf(i,:) = 1./(1j*2*pi*w - H(Sk));
        end
        
        %for i = 1:length(t)
        %    y(:,i) = 
        
        
        % integrate: x(t) = T*exp(D*t)*iT*x0 
        for i = 1:length(t)
            % We still want to drive this (linear) system - i.e.
            % x(t) = T*exp(D+(dfdu*input)*t)*iT*x0
            
            Tx = t(i)/1000;
            Ax = (T*diag(exp(d*Tx))*iT*xbar) ;
            Bu = (Tu*diag(exp(du*Tx))*iTu*in_proj(:,i));
            y(:,i) = (Ax + Bu);
            
        end
        S=[];
        
    otherwise
        % Do an actual numerical integration for a discrete epoch, if not using kernel approach
        %------------------------------------------------------------------
        for i   = 1:length(t) 
            %delays
            del = exp(P.ID).*[1 1/2 1 1 1 1 1 1];
            del = repmat(del,[1 nk]);
            del = 1./del;
            
            Tau = ( 1./(1 - (dt*del))-dt );  % <-- NOTE!
                
            if ~WithDelays 
                
                % Use a Euler integration scheme
                for j = 1:2;
                    
                    % Next step with input
                    dxdt   = f(v,drive(i,:),P,M);                

                    % delay state update: state time constants
                    dxdt = v + ( (dxdt - v)./Tau(:) );

                    % full update
                    v      = v + dt*dxdt;  

                    y(:,i) = v;
                end
                % Next step without input
                %%dxdt0   = f(v0,0,P,M);                
                %v0      = v0 + dt*dxdt0;  
                %y0(:,i) = v0;
            elseif WithDelays == 1010
                
                %delays
                del = exp(P.ID).*[1 1 1 1 1 1 1 1];
                %del = repmat(del,[1 nk]);   
                del = [del(:); ones(48,1)];
                del = 1./del;
                                
                Tau = ( 1./(1 - (dt*del))-dt );  % <-- NOTE!
                
                for j = 1:8
                    % motion
                    dxdt = f(v,drive(i,:),P,M); 

                    % delay state update
                    dx = v + (dxdt - v)./Tau(:);

                    v = v + dt*dx;
                end
                
                y(:,i) = v;
                
            elseif WithDelays == 10
                
                ff = @(x) f(x,drive(i),P,M);
                
                j = jaco(ff,spm_unvec(v,M.x),ones(size(v))/8);
                
                v = v + j\v;
                
                y(:,i) = v;
                
            elseif WithDelays == 2 
                % A simple second order integration routine akin to an
                % Euler or RK routine - but with a (matrix) delay operator
                
                for j = 1:4;%:N
                    
                    %v  = v + Q*f(spm_unvec(v,M.x),drive(i,:),P,M);
                    
                    %v0 = v0 + del'.*(Q*f(spm_unvec(v0,M.x),0.0001,P,M));                     
                    
                    % simulation integration
                    %------------------------------------------------------
                    [dx]  = f(spm_unvec(v,M.x),drive(i,:),P,M);
                    
                    %dx = Q*(dx-v);
                    %dx = Q*(dfdx*(dx-v));
                    %ddx = Q*abs(dx - v')*dt;
                    
                    %dx =  ddx*dx;
                    
                                                            
                    % augment Jacobian & delay operator, take matrix exp
                    %------------------------------------------------------
                    %e.g. expm([0   0     ]) = (expm(t*dfdx) - I)*inv(dfdx)*f
                    %          [t*f t*dfdx]
                                   
                    J = full(spm_cat({0   [];
                                    dt*dx dt*dfdx}));
                                 
                    % solve using matrix expectation & recover dv
                    %------------------------------------------------------
                    step = expm(J);
                    
                    dx = step(2:end,1);
                    
                    % delay state update
                    dx = v + (dx - v)./Tau(:);
                                        
                    v = v + dx;
                                 
                    %dx = spm_dx(dfdx,dx,diag(Q));
                    %dx = Q*dx;
                    
                    v = v + dx;
                end
                
                y(:,i) = v;%(v0 - v);
                y0(:,i) = v;%v0;
                ys(:,i) = v;
                
                
                
            elseif WithDelays == 3
                % A simple second order integration routine akin to an
                % Euler or RK routine - but with a (matrix) delay operator
                % and updating dfdx - i.e. full Jacobian integration
                
                for j = 1
                    
                    % simulation integration
                    %------------------------------------------------------
                    [dx]  = f(spm_unvec(v,M.x),drive(i,:),P,M);
                    
                    
                    % augment Jacobian & delay operator, take matrix exp
                    %------------------------------------------------------
                    J = full(spm_cat({0   [];
                                     dt*dx dt*(Q.*dfdx)}));
                                                     
                    % solve using matrix expectation
                    %------------------------------------------------------
                    step = expm(J);
                    v = v + step(2:end,1);
                    
                    % no-stim integration
                    %------------------------------------------------------
                    bx  = f(spm_unvec(v0,M.x),1e-6,P,M);
                                        
                    % augment Jacobian & delay operator, take matrix exp
                    %------------------------------------------------------
                    Jb = full(spm_cat({0   [];
                                      dt*bx dt*(Q.*dfdx)}));
                    
                    % solve using matrix expectation
                    %------------------------------------------------------
                    stepb = expm(Jb);
                    v0 = v0 + stepb(2:end,1);
                    
                end
                
                y(:,i) = v - v0; %- spm_vec(M.x);    
                
            elseif WithDelays == 44
                % Newton-Cotes quadrature, including Q
                % note this is a Lagrange polynomial method
                
                % forced temporal filter over membrane potentials -
                for j = 1:(N)
                    k1 = f(v      ,drive(i,:),P,M);
                    k2 = f(v+dt*k1,drive(i,:),P,M);
                    k3 = f(v+dt*k2,drive(i,:),P,M);
                    k4 = f(v+dt*k3,drive(i,:),P,M);
                    k5 = f(v+dt*k4,drive(i,:),P,M);

                    dxdt = ((2*dt)./45)*(7*k1 + 32*k2 + 12*k3 + 32*k4 + 7*k5);   
                    
                    %dxdt = ((2)./45)*(7*k1 + 32*k2 + 12*k3 + 32*k4 + 7*k5);   
                    
%                     J = full(spm_cat({0   [];
%                         dxdt dt*(dfdx.*Q)}));
%                     
%                     % solve using matrix expectation
%                     %------------------------------------------------------
%                     step = expm(J);
%                     v = v + step(2:end,1);
                    
                    %v = v + spm_dx(dfdx.*Q,dxdt,dt);

                    v    = v + (Q*dxdt);
                    
                    %v = v + Q*(dxdt - v);
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
                y(:,i) = v - spm_vec(M.x);  % Treat as expansion about fp          
                                
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
                y(:,i)   = (del.*Q)*v;    
            
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
            %--------------------------------------------------------------
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
%y = downsample(y,2,0);
%y0 = downsample(y0,2,0);
%y = resample(downsample(y',2,0),2,1)';

%y = lowpass(y,3,1/dt);

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
try
    yys = spm_unvec(ys,y);
catch
    yys = y;
end

series.States_without = yy;
series.States_with_inp = yys;
series.States_both = y;


% Compute the cross spectral responses from the integrated states timeseries 
%==========================================================================
[y,s,g,noise,layers] = spectral_response(P,M,y,w,npp,nk,ns,t,nf,timeseries,dt,dfdx,ci,1);

% [yn,sn,gn,noisen,layers1n] = spectral_response(P,M,yy,w,npp,nk,ns,t,nf,timeseries,dt,dfdx,ci,1);
% 
% 
% %y = exp(P.Ly(1))*abs(y - yn);
% 
% y = exp(P.Ly) * atcm.fun.moving_average(abs(y - yn),2);
% 
% series.with_inp = y;
% t = drive;
% layers = spm_unvec( spm_vec(layers1) - spm_vec(layers1n), layers1);
% %layers = layers1;


% Gaussian filtering
%---------------------------------------------------------
% if any(y)
%     [y,gm] = atcm.fun.Sig2GM(y,2);
%     y = y(:);
%     
%     % return also the gm coefficients
%     series.gm = gm;
% end
% 

end

function [y,s,g,noise,layers]=spectral_response(P,M,y,w,npp,nk,ns,t,nf,timeseries,dt,dfdx,ci,type)
% Main spectral response function with lots of options.

% Spectral Response Options
%--------------------------------------------------------------------------
DoHilbert      = 0; % take the absolute (magnitude) of the hilbert envelope
Bandpassfilter = 0; % band pass filter to [w(1)-1) : w]
DoDCT          = 0; % discrete cosine transform series before fft
IncDCS         = 0; % include semi-stochastic neuronal fluctuations       x       % ON FOR DEXPRO
DoHamming      = 0; %(1)% Cosine Hann window the model spectrum      
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
for i = 1:size(P.a,2) 
    Gu(:,i) = exp(P.a(1,i))*(w.^0);                  % P.a = constant
end

% Spectrum of channel noise (non-specific)
Gn = exp(P.b(1,i) ).*w.^(-exp(P.b(2,1))); 

% Spectrum of channel noise (specific)
for i = 1:size(P.c,2) 
    Gs(:,i) = exp(P.c(1,i) )+w.^(-exp(P.c(2,1)));     % P.c = expone
end

% Hamming to taper the edges - optimisation struggles with edge effects
%----------------------------------------------------------------------
if HamNoise
    warning off ;        % dont warn of integer operands
    H  = kaiser(nf,2.5);
    warning on;
end

% % Neuronal innovations: a discrete cosine basis set (length of P.d)
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
    
% Optional burn in/out - i.e. take transform of n:m ms instead of 0:end...
%--------------------------------------------------------------------------
burn  = atcm.fun.findthenearest(300,t); 
burno = length(t); 
if isfield(M,'burnin')
    burn = atcm.fun.findthenearest(M.burnin,t); 
end
if isfield(M,'burnout')
    burno = atcm.fun.findthenearest(M.burnout,t); 
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
    ts0       = spm_vec(J'*xseries);
        
    % add noise to make it look like a real signal
    %----------------------------------------------------------------------
    %s2n = mean(ts0(:))*1/16;
    %ts0 = ts0 + s2n*randn(size(ts0));
    %ts0 = ts0 - mean(ts0);
    
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
        warning off;
        ndmd=2;
        [Eigenvalues, Eigenvectors] = atcm.fun.dmd(yx', ndmd, dt);
        yx = Eigenvectors*yx;
        Eigenvectors = yx;
        % the user supplied j-vector no longer makes sense, ignore it -
        Ji = 1:ndmd;
        P.J(1:ndmd) = log(1.1);
        warning on;
    end
    
    warning on;
    
    % loop spatio-temporal modes (in this region) and weight them
    %----------------------------------------------------------------------
    for ij = 1:length(Ji)
        
        y0 = Eigenvectors(Ji(ij),:) ;
        Hz = w;
        
        % Window-smoothed fft
        %------------------------------------------------------------------
        warning off;
        try
            
            % State timeseries without burnin
            %--------------------------------------------------------------
            this = Eigenvectors(Ji(ij),burn:burno);
                        
            % DCT/iDCT components [for filtering only]
            %--------------------------------------------------------------
            dct_filter = 0;
            
            if dct_filter
                test = atcm.fun.adct(this)';
                pc = test;
                nn = size(test,1);
            else
                test = this;%y0;
                pc = test;
                nn = size(test,1);
            end
                        
            clear Ppf  Pfm Ppf1 Ppf2 Ppf3 

            % Select sliding-window smoothed fft or whole chunk fft
            %--------------------------------------------------------------
            UseSmooth = 0;
            
            if isfield(M,'UseSmooth') && ~isempty(M.UseSmooth)
                UseSmooth = M.UseSmooth;
            end
            
            if UseSmooth == 1
                
                % User specified FFT smoothing (num windows)
                %----------------------------------------------------------
                if isfield(M,'smth') && ~isempty(M.smth)
                    smth = M.smth;
                else
                    smth = 30;
                end

                % Smoothed Fourier transform (atcm.fun.AfftSmooth)
                %----------------------------------------------------------
                [Ppf,~,Pfm] = atcm.fun.AfftSmooth( pc, dw./dt, w, smth) ;

                nd = size(P.d,1);
                X  = spm_dctmtx(length(w),nd + 1);
                Mu = exp(X(:,2:end)*(P.d));
                
                %Ppf = atcm.fun.aenv(Ppf,6);
                
                %Ppf = Ppf(:) .* Mu(:);
                
                Ppf = idct( dct(Ppf(:)).*Mu(:) );
                
                % Retain principal eigenmode over (Gaussian) time-freq windows
                %Ppf = squeeze(Pfm)';
                %[u,s,v] = svd(Ppf');n = 1;
                %Ppf = spm_vec(u(:,n)*s(n,n)*mean(v(:,n)));
    
            elseif UseSmooth == 0 % (else use non smooth)
            
                % Non-smoothed fft version:
                %--------------------------------------------------------------
                for i = 1:nn; Ppf = atcm.fun.Afft( pc, dw./dt, w) ;end
                %Ppf = atcm.fun.HighResMeanFilt(Ppf,1,4);                
                %for i = 1:nn; Ppf = periodogram( pc(i,:), [], w,1./dt) ;end
                %Ppf = atcm.fun.HighResMeanFilt(Ppf,1,2); 
                
                %Ppf = atcm.fun.moving_average(Ppf,2);
                
                %Ppf = Ppf(:)./Gn(:);
                
                %Ppf = abs(AGenQ(Ppf))*Ppf(:);
                
                %[~,GL] = AGenQ(Ppf);
                
                %Ppf = GL*Ppf(:);
                
                nd = size(P.d,1);
                X  = spm_dctmtx(length(w),nd + 1);
                Mu = exp(X(:,2:end)*(P.d));
                
                %Ppf = atcm.fun.aenv(Ppf,6);
                
                Ppf = Ppf(:) .* Mu(:);
                
                %Ppf = idct( dct(Ppf(:)).*Mu(:) );
                
            elseif UseSmooth == 2

                K = 24;
                                
                % Compute Dynamic Mode Decomposition on TF data
                [Ppf,~,c] = atcm.fun.tfdecomp(pc,dt,w,4,1,@median);
                
                %c = AGenQ(Ppf)*c;
                                
                % Remove ~1/f system noise
                %c = abs(c) ./ repmat(Gn(:),[1 size(c,2)]);
                                                
                % Eigenvectors of the DMD of TF matrix are spectra of modes
                [Eigenvalues, ev, ModeAmplitudes, ModeFrequencies, GrowthRates, POD_Mode_Energies] = atcm.fun.dmd(c, K, dt);
                
                %Pp  = c*ev';
                %Ppf = spm_vec(max(abs(Pp)'));
                
                % Compute modal frequencies
                ModeFrequencies=(angle(Eigenvalues)/pi)*(1/(2*dt));
                                
                % Make values positive (as per a psd estimate)
                ModeAmplitudes  = abs(ModeAmplitudes);
                ModeFrequencies = abs(ModeFrequencies);
                ModeFrequencies = round(ModeFrequencies*100)./100;
                
                %[ModeFrequencies ModeAmplitudes]
                
                % remove duplicates
                [~,I]=unique(ModeFrequencies);
                
                ModeFrequencies = ModeFrequencies(I);
                ModeAmplitudes  = ModeAmplitudes(I);
                
                % Remove stuff outside freqs of interest
                G = find( ModeFrequencies>w(1) & ModeFrequencies<w(end) );
                
                %wt = round(ModeFrequencies)/max(w);
                %ModeAmplitudes = ModeAmplitudes(:).*wt(:);
                                                               
                % Convert to (Gaussian) series and average
                for ijf = 1:length(G)
                   PF(ijf,:) = atcm.fun.makef(w,ModeFrequencies(G(ijf)),ModeAmplitudes(G(ijf)),1.5);
                end
                 
                 % average and smooth
                 Ppf = spm_vec(mean(PF));
                 wt  = w(:)./(max(w)*exp(P.Ly));
                 Ppf = Ppf.*( wt );
                
                %[~,GL] = AGenQ(Ppf);                
                %Ppf = GL*Ppf(:);
                
                % dct transform
                %nd  = size(P.d,1);
                %X   = spm_dctmtx(length(w),nd + 1);
                %Mu  = exp(X(:,2:end)*P.d);
                %Ppf = idct( dct(Ppf(:)).*Mu(:) );
                                
            elseif UseSmooth == 3
                
                [Ppf,~,c] = atcm.fun.tfdecomp(pc,dt,w,8,   1,@max);
                
                [~,GL] = AGenQ(std(c'));
                
                Ppf = GL*Ppf(:);
                %Ppf = AGenQ(Ppf)*Ppf;
                %Ppf = AGenQ(std(c'))*Ppf;
                %Ppf = AGenQ(Ppf)*Ppf;
                
                nd  = size(P.d,1);
                X   = spm_dctmtx(length(w),nd + 1);
                Mu  = exp(X(:,2:end)*P.d);
                Ppf = idct( dct(Ppf(:)).*Mu(:) );
                
            elseif UseSmooth == 4
                 
                Ppf = atcm.fun.Afft( pc, dw./dt, w)' ;
                Ppf = atcm.fun.moving_average(Ppf,2);
                
                [~,~,c] = atcm.fun.tfdecomp(pc,dt,w,12,2,@mean);
                
                C = cov(c');
                                
                % Laplacian smoothing
                A  = real(C) .* ~eye(length(C));
                N  = size(A,1);
                GL = speye(N,N) + (A - spdiags(sum(A,2),0,N,N))/4;
                
                Ppf = GL*Ppf;
                Ppf = GL*Ppf;
                Ppf(Ppf<0) = 0;
                %Ppf = abs(Ppf);
                
                                
                nd  = size(P.d,1);
                X   = spm_dctmtx(length(w),nd + 1);
                Mu  = exp(X(:,2:end)*P.d);
                Ppf = idct( dct(Ppf(:)).*Mu(:) );
                
                
            end
            
            % De-NaN/inf the spectrum
            %--------------------------------------------------------------
            Pf = Ppf;            
            Pf(isnan(Pf))=0;
            Pf(isinf(Pf))=0;
            
            % return the series / orthogonal VMD components if selctd
            %--------------------------------------------------------------
            layers.ssa_pc{ins,ij,:,:} = pc;
            layers.pst_burn = t(burn:end);
            
            % Allow a gain specifically for imaginary components (for CSDs)
            %--------------------------------------------------------------
            if isfield(P,'iL');
                Pf = real(Pf) + sqrt(-1)*exp(P.iL(ins))*imag(Pf);
            end
            
            % make sure its a nx1 vector
            %--------------------------------------------------------------
            Pf  = spm_vec(Pf);
            %Pf(Pf<0) = -Pf(Pf<0);
                        
        end
    
        % Make it non-sparse and vectorised
        warning on;
        try
            Pf = (Pf)';
        catch
            Pf = inf*size(w);
            y = Pf;
            s = timeseries;
            g = ts;
            return;
        end
        
        Pf = full(Pf)';
        J  = full(exp(P.J));
        
        if DoHamming
            H = .5+hamming(nf,'periodic');
            %[val,ind] = max(H);
            %H(1:ind)  = val;
            %H = H - min(H);
            Pf = Pf(:).*H(:);
        end

        % store the weighted and unweighted population outputs
        %------------------------------------------------------------------
        layers.unweighted(ins,ij,:) = ( Pf             )     ;% * exp(real(P.L(ins)));
        layers.weighted  (ins,ij,:) = ( Pf * abs(J(Ji(ij))) );% * exp(real(P.L(ins)));
        layers.iweighted (ins,ij,:) = ( Pf * abs(J(Ji(ij))) );% * exp(real(P.L(ins)));
        
    end   % end of state contribution(s) to region loop
end   % end of loop of regions


clear Pf

% Now compute node proper CSDs from sum of weighted cells [modes] per region
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
%Pf = (Pf)/length(w);


% Addition of system noise & Lead field scaling: L, Gu, Gn, Gs
%--------------------------------------------------------------------------
for ins = 1:ns
        
    % Peak detection and multivariate distribution fitting
    %----------------------------------------------------------------------
%     if isfield(M,'threshold') && ~isempty(M.threshold)
%         threshold = M.threshold;
%     else
%         threshold = 0.8;
%     end
%         
%     %estimate system noise (exponential decay) & remove before fit
%     Pf(:,ins,ins) = atcm.fun.HighResMeanFilt(Pf(:,ins,ins),1,2);
%     
%     modelexp   =  'exp(a)*(x.^(-exp(b)))';
%     aperiod    = atcm.fun.c_oof(w,(Pf(:,ins,ins)),modelexp);
%     
%     
%     %Pf(:,ins,ins) = max(real([Pf(:) aperiod(:)]'));
%     
%     %dfn        = real(Pf(:,ins,ins))./sum(real(Pf(:,ins,ins)));
%     dfn        = (Pf(:,ins,ins));
%     dfn        = dfn(:) - aperiod(:);
%     sdfn       = sum(dfn);
%     dfn        = dfn./sum(dfn);
%     
%     %dfn  = detrend(dfn);
%         
%     [odfn,ord] = sort(dfn,'descend');
%         
%     n    = atcm.fun.findthenearest(cumsum(odfn),threshold); 
%     pk   = ord(1:n);
%      wint = 1:length(w);
%         
%     ResY = squeeze(M.y{ci}(:,ins,ins));
    
    % Fitting, smoothing and/or filtering of the LFP spectrum for region i
    %----------------------------------------------------------------------
    if UseSmooth == 0
        wt   = rescale(w(:).^-1);
    else
        wt = w'.^0;
    end
    wt = w'.^0;
    
    if isfield(M,'dist_method') && ~isempty(M.dist_method)
        dist_method = M.dist_method;
    else
        dist_method = 4;
    end
    
    switch dist_method
        case 1
            % smoothed spline
            Pf0   = csaps(w(pk),Pf(pk,ins,ins).*wt(1:n),0.01,w);
        case 2
            % normal cubic spline
            Pf0   = spline(w(pk),Pf(pk,ins,ins).*wt(1:n),w);
        case 3
            % gaussians
            
            % Identify peaks using diff(x)
            Pfx = Pf(:,ins,ins);
            d = diff(Pfx);
            for i = 1:length(Pfx)-1; 
                if d(i) > 0 && d(i+1) < 0; 
                     pk(i+1) = i; 
                else pk(i+1) = 0;
                end;
            end
            
            pk = find(pk);
                        
            % generate initial gaussians
            for i = 1:length(pk)
               Pf0(i,:)   = atcm.fun.makef(wint,wint(pk(i))-1,Pf(pk(i),ins,ins),2.6,'gaussian');
            end
            
            % generate a width (component proportion) function
            wfun  = @(x) atcm.fun.makef(wint,wint(pk)-1,Pf(pk,ins,ins),x,'gaussian');
            gfun  = @(x) sum( (Pfx - spm_vec(wfun(x))).^2 );
            initw = ones(size(pk));
            
            % solve with fminsearch
            [X,F] = fminsearch(gfun,initw);
            
            Pf0 = wfun(X);
            
            %Pf0 = max(Pf0);
            %Pf0 = atcm.fun.moving_average(Pf0,2);
            
            %[u,s,v] = spm_svd(Pf0);
            %Pf0 = sum( abs(u(:,1:4)'*Pf0) ,1);
            
        case 4
            % optimised gaussian, cauchy, laplacian or gamma
            Pf0   = atcm.fun.findbestdist(wint,wint(pk)-1,Pf(pk,ins,ins).*wt(1:n),2.6*ones(length(find(pk)),1),ResY);
            %Pf0   = atcm.fun.findbestdist(wint,wint(pk)-1,Pf(pk,ins,ins).*wt(1:n),ones(length(find(pk)),1),ResY);
        case 5
            % linear interpolation with smoothing
            Pf0 = interp1([w(1)-1 w(pk) w(end)+1],Pf([1 pk' end],ins,ins).*[1; wt(1:n); 1],w,'pchip');
            %Pf0 = atcm.fun.HighResMeanFilt(Pf0,1,4);
        case 6
            %Pf0 = atcm.fun.HighResMeanFilt(Pf(:,ins,ins),1,2);
            %Pf0 = atcm.fun.tsmovavg(Pf(:,ins,ins)','e',4);
            %Pf0 = atcm.fun.moving_average(Pf(:,ins,ins),2);
            Pf0 = Pf(:,ins,ins);
            %Pf0 = sdfn*dfn;
            %Pf0 = atcm.fun.moving_average(Pf0,2);
            %Pf0 = atcm.fun.HighResMeanFilt(Pf0,1,2);
            
%             k = 8;
%             warning off;
%             try
%                 theta = fitgmdist([w(:) real(Pf0)],k);
%                 Pf0 = makef(w(:),theta.mu(1:k),theta.mu(k+1:end),squeeze(theta.Sigma(1,1,:)));
%             catch
%                 Pf0 = ones(size(Pf0))*1e20;
%             end
%             warning on;
%             Pf0 = makef(w(:),theta.mu(1:k),theta.mu(k+1:end),squeeze(theta.Sigma(1,1,:)));
            
        case 7
            Pf0 = lowpass(squeeze(Pf(:,ins,ins)),.2,1);
        case 8
            Pf0 = squeeze(Pf(:,ins,ins));
            curv = atcm.fun.c_oof(w,Pf0);
            Pf0 = Pf0 - curv;
            m   = fit(w.',real(Pf0),'Gauss5');
            Pf0 = curv + m(w);
        case 9
            Pf0 = atcm.fun.component_spectrum(w,squeeze(Pf(:,ins,ins)),5);
 
    end
        
%     cv = atcm.fun.estcov(Pf0(:),length(Pf0));
%     [u,d] = eig(cv);
%     d = diag(d);
%     [~,I]=sort(d,'descend');
%     u = u(:,I);
%     d = d(I);
%     Pf0 = u(:,1)'*cv;
    
    
    % Add aperiodic component back in
    %Pf(:,ins,ins) = aperiod(:) + Pf0(:) ;
    Pf(:,ins,ins) = Pf0(:) ;
    
    % Electrode gain & smoothing
    %----------------------------------------------------------------------
    Pf(:,ins,ins) = exp(P.L(ins))*Pf(:,ins,ins);

    % Add noise to (in frequency space) to this LFP channel spectrum
    %----------------------------------------------------------------------
    addnoise=0;
    if addnoise
        % Multiply in the {dct} neuronal fluctuations - Gu [P.a / P.d]
        for i = 1:length(Hz)
            Pf(i,ins,ins) = (sq(real(Pf(i,ins,ins)))*real(diag(Gu(i,ins)))*sq(real(Pf(i,ins,ins)))) + sqrt(-1)*imag(Pf(i,ins,ins));
        end
    end

end     

% Recompute cross spectrum of regions (CSDs) incl. noise
%--------------------------------------------------------------------------
for inx = 1:ns
    for iny = 1:ns
        if inx ~= iny
            Pf(:,inx,iny) = squeeze( Pf(:,inx,inx) ) .* conj( Pf(:,iny,iny) ) ;
        end
    end
end

% Re-incorporate other noise components for auto (Gs) and cross (Gn) spectra
%--------------------------------------------------------------------------
addnoise=0;
if addnoise        
    for i = 1:ns
        for j = 1:ns
            % Autospectral noise / innovations [P.b]
            Pf(:,i,j) = (Pf(:,i,j).*Gs(:,i)) + Gn(:);
            if j ~= i
                % Cross spectral noise / innovations [P.c]
                %Pf(:,j,i) = Pf(:,j,i) .* Gs(:,i);
                %Pf(:,i,j) = Pf(:,j,i);
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
if ns == 1 
    y = real(Pf);
    %y=Pf;
    %y=abs(Pf);
else
    y = Pf;
end

s = timeseries;
g = ts;

end

% Helper subfunctions

function g = thresh_fft(x,Ppf,pk,M,ci,ins,w)
pk  = Ppf>(mean(Ppf)+std(Ppf)./x);
Ppf = atcm.fun.makef(w,w(pk)-1,Ppf(pk),ones(length(find(pk)),1)*2);
g = sum( (spm_vec(squeeze(M.y{ci}(:,ins,ins))) - Ppf(:) ).^2 );
end

function [x] = sq(x)
% squeeze function for csds
if size(x,3) > 1, x = squeeze(x); else, x = x(:); end

end

function j = adfdx(IS,P,M,order,xhat)
% faster computation of the jacobian (or hessian) states matrix
% AS
warning off ;

if nargin < 5 || isempty(xhat);
    xhat = spm_vec(M.x);
end

fx    = spm_vec(feval(IS,xhat,0,P,M));
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



