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
% [y,w,s,g,t,pst,layers,noise,firing,QD,Spike] = atcm.integrate_1(P,M,U)
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
% or Runge-Kutta scheme to first generate a states time series 
% (voltages, currents). Finally, weighting & an fft of  the resultant 
% timeseries provide the spectral response. This is nice because one can access 
% % both the continuous and discrete time responses of the system. 
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
% By default a modified Euler scheme incorporating a delay operator is used. 
% To switch from a numerical integration to another option:
% Set M.IntMethod = 'kernels'   ... use spm kernel scheme
%       ~~~~      = 'ode45'     ... built in matlab ode solver
%                 = ' '         ... use a numerical method (default)
%                 = 'linearise' ... local linearisation using
%                                   eigendecomposition of the jacobian
% 
% If using a numerical method, switch the method type:
% Set M.intmethod = 0  ... Euler no delays
%       ~~~~      = 2  ... Euler with delays
%                 = 44 ... Newton-Cotes with delays
%                 = 21 ... integration with bilinear jacobian
%                 = 23 ... full jacobian integration
%                 = 24 ... stochastic equation integration
%                 = 8  ... 8th-order Runge-Kutta w/ delays
%                 = 45 ... 4th order Runge-Kutta w/ delays *
%
% The RK4 DDE implementation is as follows:
% 
%  k1 = Q*f(v          ,drive(i)     ,P,M);
%  k2 = Q*f(v+0.5*dt*k1,drive(i)+dt/2,P,M);
%  k3 = Q*f(v+0.5*dt*k2,drive(i)+dt/2,P,M);
%  k4 = Q*f(v+    dt*k3,drive(i)     ,P,M);
% 
%  dxdt = (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
%  v         = v + dxdt;
%
% Also required: SPM12 w/ DCM, plus aoptim/AO.m for param optimimsation.
% Dr Alexander D Shaw | 2020 | alexandershaw4[@]gmail.com


% w, initial states, dt | fs is specified & generate sampletimes & inputs
%--------------------------------------------------------------------------
w     = M.Hz;                     % FoI (w)
x     = M.x;                      % model (hidden) states
Kx    = x;                        % pre-fp x, for kernel scheme
dt    = 1/600;                    % sample rate
Fs    = 1/dt;                     % sampling frequency
tn    = 3;                        % sample window length, in seconds
pst   = 1000*((0:dt:tn-dt)');     % peristim times we'll sample

% stop P from becoming complex
%P = spm_unvec( real(spm_vec(P)), P);

P.L = real(P.L);
P.J = real(P.J);

if nargin == 4 && ~isempty(varargin{1})
    %fprintf('using user-supplied M.x\m');
    M.x = spm_unvec(varargin{1},M.x);
    x = M.x;
end

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

        drive = drive(:)';%.*NS(:); 
    case 1
        
        % create a noise distribution in frequency space for endogenous inp
        fdrive =  1./(1j*2*pi*w - (50)*exp(P.R(1)));
        drive = ifft(fdrive,length(pst)*2);
        drive = exp(P.R(2))*drive(1:length(pst))*5000;

        % and an ERP like input for exogenous input to thalamic
        %fre = 40 * exp(P.R(2));
        %drive = exp(P.R(1)) * sin(2*pi*fre*pst./1000);
                  
    case 2
        % For ERP inputs...
        %------------------------------------------------------------------
        delay  = 6 * exp(P.R(1));             % bump
        scale1 = 8  * exp(P.R(2));
        drive  = atcm.fun.makef(pst,delay,scale1,16);

        offset = max(drive)/8;
        intercept = atcm.fun.findthenearest(offset,drive);
        drive(intercept(1):end) = offset;

        drive=drive';
        
    case 3
        a = 1 * exp(P.R(1));
        b = 4 * exp(P.R(2));
        drive =  a + (b-a).*rand(1,length(pst));        
    case 4
        % TWO oscillatory inputs...
        %------------------------------------------------------------------
        mu1   = 0.5*exp(P.R(1));                % mean amplitude
        mu2   = 1.0*exp(P.R(2));
        mf1   = 100*exp(P.R(3));                  % frequency
        mf2   = 1*exp(P.R(4));
        
        drive(:,2) = (mu1 * sin(2*pi*mf1*(pst/1000)) ) ;  
        drive(:,1) = (mu2 * sin(2*pi*mf2*(pst/1000)) );

        drive = prod(drive,2);
        
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
        
        % For oscillatory inputs...
        %------------------------------------------------------------------
        mu    = .2*exp(P.R(1));                      % mean amplitude
        mf    = 10*exp(P.R(2));                      % frequency
        drive(:,2) = mu * ( (sin(2*pi*mf*(pst/1000))) );%...

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

    case 7

        drive = zeros(56,length(pst));
        If = [10 55 80 40 50 20 10 10];
        for i = 1:8
            drive(i,:) = exp(P.R(i)) * ( (sin(2*pi*(If(i))*(pst/1000))) );
        end

        % break phase locking between layers
        d = [1 2 2 3 3 4 0 0];
        for i = 1:8
            ni = round((d(i)/200)/dt);
            drive(i,:) = [zeros(1,ni), drive(i,1:end-ni)];
        end

        drive = sum(drive,1);
    
    case 8
         W = atcm.fun.hnoisebasis(length(w),exp(P.R(1:2)));
         [spec] = atcm.fun.asinespectrum(w,pst);
         drive = (W')*shiftphase(spec,6);
    case 9

        %  retinogeniculate oscillations are broadband incl. gamma;
        % https://www.frontiersin.org/articles/10.3389/neuro.06.004.2009/full
        x0 = pst/1000;
        w0  = exp(P.d(1));
        a0 = exp(P.d(2));
        a1 = exp(P.d(3));
        
        b1 = exp(P.d(4));

        a2 = exp(P.d(5));
        b2 = exp(P.d(6));

        drive =  a0 + a1*cos(x0*w0) + b1*sin(x0*w0) + ...
              a2*cos(2*x0*w0) + b2*sin(2*x0*w0);

        drive = drive';


end

if isfield(M,'input')
    dl = length(drive);
    drive = M.input(1:dl);
end

fso=[];

if ~isfield(M,'timefreq')
    M.timefreq = 0;
end

% expansion (fixed) point: trial & parameter effects are deviations from here
%--------------------------------------------------------------------------
f    = spm_funcheck(M.f); 

% solve for a fixed point, or not
if solvefp; 
    %x    = atcm.fun.solvefixedpoint(P,M);
        %x = atcm.fun.solvefixedpoint(P,M,[],-70);
        x = spm_unvec(atcm.fun.alexfixed(P,M,1e-6,[],1),M.x);
        
else ;      x    = x;
end

%DM = struct; DM.M=M;DM.pE=P;
%x = FindSteadyState(DM,128,1/1200);

M.x  = x;
v    = spm_vec(x);
NoFX = 0;

% flag no modulatory effects in this model
if isempty(U) || (isfield(U,'X') && isempty(U.X))
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
    [y{c},w,s{c},g{c},t{c},layers{c},firing{c},QD{c},Spike{c},condel{c},series{c}] = ...
        dodxdt(pst,f,v,Q,M,dt,w,drive,Kx,U,method,solvefp,c,fso);

end

% outputs
%--------------------------------------------------------------------------
other.firing = firing;
other.QD = QD;
other.Spike = Spike;
other.drive = drive;
other.condel=condel;
other.series = series;
other.dt = dt;
other.Fs = Fs;

end

function [y,w,s,g,t,layers,firing,QD,spike,condel,series] = ...
                            dodxdt(t,f,v,P,M,dt,w,drive,Kx,U,method,solvefp,ci,fso)
% Numerical integration, signal weighting and FFT of timeseries

% Choose either the DCM kernels scheme, or another dx = f(x) update scheme
%--------------------------------------------------------------------------
if isempty(method)
    method = 45; % default to the delayed RK45
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
[fx, dfdx] = f(M.x,0,P,M);
condel = 1;
QD     = 1;

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

% setup a baseline jacobian using finite differecnes
fun = @(v) f(spm_unvec(v,M.x),drive(1),P,M);
Jf = jaco(fun,v,ones(size(v))/128,0,2);
vi = v*0;

xb = M.x(:);

switch IntMethod

    case 'amtf'

        % another transfer function based on the eigenspectrum of the
        % numerical jacobian
        y = exp(P.J)*amtf(P,M,U);

    case 'spm';
        x;
        P.J=P.J';
        M.dt = dt;
        [~,y] = spm_int_L(P,M,drive);
        S=[];

    case 'kernels'
        
        % Just use the kernels approach?
        %------------------------------------------------------------------
        [K0,K1,K2,H1] = spm_kernels(M,P,length(t),dt);
        y = H1(:,1:end-1)';
        S = [];
        
    case 'ode45'

        [~,~,D] = f(v,drive(1),P,M);
        Q = eye(length(D)) - (D);

        % matlab build in ode solver
        %ode = @(t,v,P,M,f) spm_vec( Q*f(spm_unvec(v,M.x),t,P,M) );
        %ode = @(t,v,P,M,f) spm_vec( Q*f(spm_unvec(v,M.x),drive(1),P,M) );
        
        ode = @(t,v) spm_vec( Q*f(spm_unvec(v,M.x),t,P,M) );

        opts = odeset;
        %opts.MaxStep = 1/600;
        
        [~,y]   = ode113(ode,t/1000,spm_vec(v));%,drive(1),P,M,f);
                
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
        dfdx    = dfdx./norm(full(dfdx));
                            
        % compute df/du (inputs)
        %for i = 1;%:length(drive)
            dfdu    = spm_diff(M.f,M.x,drive(1)*dt,P,M,2);  
        %end
        %dfdut = dfdu;
        %dfdu  = sqrt(dfdu*dfdu');
        %dfdu  = dfdu./norm(full(dfdu));
        
        % eigenvectors and values - dfdx
        [T,D]   = eig(full(real((dfdx))));
        iT      = pinv(T);
        d       = diag(D);
        
        dgdv  = T;
        dvdu  = pinv(T).*dfdu;
                        
        % integrate: x(t) = T*exp(D*t)*iT*x0 
        for i = 1:length(t)            
            Tx = t(i)/1000;
            Ax = (T*diag(exp(d*Tx))*iT*xbar) ;
            
            if i == 1
                v0 = dt*Ax;
            end

            v0 = v0 + dt*Ax;
            y(:,i) = v0;%(Ax);% + Bu);
            
        end
        S=[];
        
    otherwise
        % Do an actual numerical integration for a discrete epoch, if not using kernel approach
        %------------------------------------------------------------------
        for i   = 1:length(t) 
            
            if ~WithDelays 
                                
                % Use a Euler integration scheme
                for j = 1;
                    
                    % Next step with input
                    [dxdt] = f(v,drive(i,:),P,M);      

                    % full update
                    v      = v + dt*dxdt;
                    y(:,i) = v;
                end
                
            elseif WithDelays == 2
                
                % TWO-point RK method
                k1  = f(v,drive(i,:),P,M);  
                k2  = f(v+dt*k1,drive(i,:),P,M);
                phi = 0.5*k1 + 0.5*k2;              
                v = v + dt*phi;

                % State Delays - interpolated
                %--------------------------------------------------
                d = 100*[.006 .002 .001 .004 .001 .008 .001 .008].*exp(P.ID);
                d = repmat(d,[1 nk]);
                L = (d);

                for j = 1:length(L)
                    ti = real(L(j))/dt;
                    if i > 1 && any(ti)
                        pt = t(i) - ti;
                        if pt > 0
                            v(j) = interp1(t(1:i), [y(j,1:i-1) v(j)]', pt);
                        end
                    end
                end

                y(:,i) = v;
                
                
            elseif WithDelays == 1
                 
                % Matrix delays and rate (integration) constants
                %--------------------------------------------------
                if i == 1
                    [k1,J,D] = f(v,0*drive(:,i),P,M);
                end

                D   = real(D);
                Q   = (1 - D).*(~~real(J));
                QJ  = Q.*J;
                ddt = dt;

                if i > 1
                    % endogenous noise via AMPA receptors
                    %--------------------------------------------------
                    d       = (drive(i) - drive(i-1))./exp(P.a(1));
                    v(9:16) = v(9:16) + d;
                end

                R=P;

                % 
                %--------------------------------------------------
                dxdt  = dt*f(v,0,P,M);
                dxdt1 = dt*f(v + dt*dxdt,0,P,M); 
                g    = (dxdt + dxdt1)./2 ;
                b    = J\g;
                dxdt = QJ*b ;
                %dxdt = spm_dx(dt*QJ,g);

                v    = v + dxdt;

                % Full update
                %--------------------------------------------------
                y(:,i) =   (v - xb);

                     
                
            elseif WithDelays == 1.5
                
                % *backward* euler
                if i == 1
                    [dxdt] = f(v,drive(i,:),P,M);    
                    v      = v + dt*dxdt;
                    y(:,i) = v;
                else
                    % Next step with input
                    [dxdt] = f(v,drive(i,:),P,M);    
                    %v = v + dt*dxdt;                    
                    y(:,i) = y(:,i-1) + dt*dxdt;
                    v = y(:,i);
                    
                end
                               
            elseif WithDelays == 1010
                
                %delays
                del = exp(P.d(:)').*[1 1 1 1 1 1 1 1];
                del = repmat(del,[1 nk]);
                del = 1./del;
                Tau = ( 1./(1 - (dt*del))-dt );  % <-- NOTE!
                dxdt = f(v,drive(i,:),P,M);

                % delay state update
                dx = v + (dxdt - v)./Tau(:);
                v = v + dt*dx;

                y(:,i) = v;
                
            elseif WithDelays == 10
                
                ff = @(x) f(x,drive(i),P,M);  
                j = jaco(ff,spm_unvec(v,M.x),ones(size(v))*exp(-8));
                
                warning off;
                v = v + j\v;
                warning on;
                
                y(:,i) = v;
                
            elseif WithDelays == 33
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
                
                k1 = f(v      ,drive(i,:),P,M);
                k2 = f(v+dt*k1,drive(i,:),P,M);
                k3 = f(v+dt*k2,drive(i,:),P,M);
                k4 = f(v+dt*k3,drive(i,:),P,M);
                k5 = f(v+dt*k4,drive(i,:),P,M);
                
                dxdt = ((2*dt)./45)*(7*k1 + 32*k2 + 12*k3 + 32*k4 + 7*k5);
                v    = v + (dxdt);
                
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
                    
                    % RK45 with delays
                    if i > 1
                        
                        % Matrix delays and rate (integration) constants
                        %--------------------------------------------------
                        if i == 2
                            D   = real(D);
                            Q   = (1 - D*dt);%.*(~~real(J));%inv(eye(npp*nk) - D);
                            QJ  = Q.*J;
                            ddt = dt;
                        end  

                        % Function that implements delays
                        delf = @(dxdt,v) v + QJ*(J \ (dxdt - v) );

                        % endogenous inputs
                        %--------------------------------------------------
                        v(16) = v(16) + exp(P.a(1));
                        v(12) = v(12) + exp(P.a(2));
                        v(10) = v(10) + exp(P.a(3));
                        v(18) = v(18) + exp(P.a(4));
                        v(19) = v(19) + exp(P.a(5));
                        v(26) = v(26) + exp(P.a(6));

                        R=P;
                        
                        drive(i) = 0;

                        % integrate w 4-th order Runge-Kutta method.
                        %--------------------------------------------------
                        k1 = f(v             ,drive(i),R,M);

                        k2 = f(v+0.5.*ddt.*k1,drive(i),R,M);

                        k3 = f(v+0.5.*ddt.*k2,drive(i),R,M);

                        k4 = f(v+     ddt.*k3,drive(i),R,M);
                        
                        dxdt = (ddt/6).*(k1 + 2*k2 + 2*k3 + k4);

                        if i > 2
                            % from the update dx = f(x), can can recover
                            % which x lead to which change in dx assuming a
                            % static Jacobian map; 
                            %     dv = x + dx
                            %     b  = J \ dx
                            %     dv = x + J*b
                            %     dv = x + (Q*J)*b  <-- delays in state flow

                            g    = dxdt - v;
                            b    = J\g;
                            dxdt = v + QJ*b ;

                            dxdt = dt*f(v+dxdt,drive(i),R,M);
                        end
                        
                        v    = v + dxdt;
                                               
                        % Full update
                        %--------------------------------------------------
                        y(:,i) =   (v );

                    else

                        % 4-th order Runge-Kutta method
                        %--------------------------------------------------
                        [k1,J,D] = f(v    ,0*drive(:,i),P,M);
                        k2 = f(v+0.5*dt*k1,0*drive(:,i),P,M);
                        k3 = f(v-0.5*dt*k2,0*drive(:,i),P,M);
                        k4 = f(v+    dt*k3,0*drive(:,i),P,M);

                        dxdt      = (dt/6)*(k1+2*k2+2*k3+k4);
                        v         = v + dxdt;
                        y(:,i)    = v ;%- xb;

                    end


            end  
            
            % firing function at dxdt - the sigmoid
            %--------------------------------------------------------------
            VR  = -52;
            V   = spm_unvec(v,M.x);            
            
            R  = 2/3 * exp(P.S);
            FF = 1./(1 + exp(-R'.*(v(1:8)-VR)));

            RS = 30 ;
            Fu = find( v(1:8) >= VR ); FF(Fu) = 1;
            Fl = find( v(1:8) >= RS ); FF(Fl) = 0;
            
            m(i,:)  = FF;

            
            % log whether membrane potential crossed threshold
            %--------------------------------------------------------------
            try
                fired     = find(squeeze(V(:,:,1)) >= VR);
                firings   = [firings; [i+0*fired',fired'] ];
            catch
                % doesn't work for multi-node models
                fired=[];
                firings=[];
            end
        end
end

warning on;

S = m;


% Reshape to model state space outputs
%--------------------------------------------------------------------------
[ns,npp,nk] = size(M.x);
y(isnan(y)) = 0;
y(isinf(y)) = 0;

% recompute dfdx now we've propoerly burned in 
[fx, dfdx] = f(spm_unvec(y(:,end),M.x),0,P,M);

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
[y,s,g,layers] = spectral_response(P,M,y,w,npp,nk,ns,t,nf,timeseries,dt,dfdx,ci,1,fso,drive);

%PP = P;
%PP.J = [1 1 1 1 1 0 0 0];
%firy = spectral_response(PP,M,m',w,npp,nk,ns,t,nf,timeseries,dt,dfdx,ci,1,fso,drive);

% jf = @(x)struct('type','.','subs',J)
% jd = @(x) spectral_response(jf(x),M,y,w,npp,nk,ns,t,nf,timeseries,dt,dfdx,ci,1,fso,drive)

%y = y + firy;

end

function [y,s,g,layers]=spectral_response(P,M,y,w,npp,nk,ns,t,nf,timeseries,dt,dfdx,ci,type,fso,drive)
% Main spectral response function with lots of options.

% Spectral Response Options
%--------------------------------------------------------------------------
DoHamming      = 0; %(1)% Cosine Hann window the model spectrum      

if isfield(M,'DoHamming')
    DoHamming = M.DoHamming;
end

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
% for ins = 1:ns
% 
%     % J-weighted sum of this channel / region
%     %----------------------------------------------------------------------
%     xseries   = full( squeeze(y(ins,:,:,:)) );
%     xseries   = reshape(xseries, [npp*nk,length(t)] );
%     ts0       = spm_vec(J'*xseries);
% 
%     % apply electrode gain & store this channel / node
%     %----------------------------------------------------------------------
%     ts(ins,:) = ts(ins,:) + ts0' ;
% 
% end
ts=0;


% adjustment for when frequency intervals ~= 1 Hz; i.e. dt = dw/dt
%--------------------------------------------------------------------------
dw = 1./(w(2)-w(1));

% This is the main loop over regions to calculate the spectrum
%==========================================================================
for ins = 1:ns
        
    % extract time series of all states from this region
    %----------------------------------------------------------------------
    if ndims(y) == 4
        yx = reshape( squeeze(y(ins,:,:,:)), [npp*nk,length(t)]);
    else
        yx = y;
    end
    
    yx = yx(:,burn:burno);

    % use eigenspectrum of numerical Jacobian
    %----------------------------------------------------------------------
    UseSS = 0;
    if UseSS
        C = cov(yx');
        D = -C.*dfdx'/2;
        [v,s]=eig(full(D),'nobalance');
        s=diag(s);
        Kw = spm_s2csd(s,w);
        Ppf = sum(Kw(:,1:8),2);
        P.J=1;
        layers.iweighted(1,1,:) = Ppf;
        continue;
    end

    if isfield(M,'dmd') && M.dmd
        [evec,eval] = atcm.fun.dmd(yx', 1, dt);
        yx = eval*yx;
        P.J = 1;
    end


    % loop spatio-temporal modes (in this region) and weight them
    %----------------------------------------------------------------------
    for ij = 1:length(Ji)
        
        y0 = yx(Ji(ij),:) ;
        Hz = w;
        
        % Spectral responses of (observable) states
        %------------------------------------------------------------------
        try
            
            clear Ppf  Pfm Ppf1 Ppf2 Ppf3    

            ys  = yx(Ji(ij),:);
            ys  = atcm.fun.bandpassfilter(ys,1/dt,[w(1) w(end)]);
            
         
            if ij == 1
                F   = dftmtx(size(ys,2));   
                N   = length(F);
                f   = (1/dt) * (0:(N/2))/N;
                Mt  = pinv(atcm.fun.cdist(f',w'))';

                if isfield(M,'GFFTM') && ~isempty(M.GFFTM)
                    Mt = M.GFFTM; 
                else
                    % generate the interpolation matrix s.t. Gauss Kernel
                    % smoothing on both axes
                    %GM = atcm.fun.VtoGauss(ones(length(w),1),1);
                    %FM = atcm.fun.VtoGauss(ones(length(f),1),1);
                    for i = 1:size(Mt,1)
                        %Mt(:,i) = atcm.fun.VtoGauss(Mt(:,i),8)*Mt(:,i);
                        Mt(i,:) = abs(Mt(i,:))*atcm.fun.VtoGauss(abs(Mt(i,:)));
                    end
                    %Mt = FM*Mt;
                end
            end          

            %for i = 1:length(w); I(i) = atcm.fun.findthenearest(f,w(i)); end
            
            % Pf(w) = fft(y)*iM, where iM is the inverse distance between
            % natural frequency vector and those of interest, subject to
            % each column of iM conforming to a Gaussian
            Ppf = (ys*F(:,1:length(f))*Mt)./N;                  
            Ppf = abs(Ppf);

            Ppf = atcm.fun.agauss_smooth(Ppf,1); %  kernel smoothing
            %Ppf = gau_signal_decomp(Ppf,4);

           
            % De-NaN/inf the spectrum
            %--------------------------------------------------------------
            Pf = Ppf;            
            Pf(isnan(Pf))=0;
            Pf(isinf(Pf))=0;
                                   
            % make sure its a nx1 vector
            %--------------------------------------------------------------
            Pf  = spm_vec(Pf);
                        
        end
    
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
            H  = .5+hamming(nf,'periodic');
            Pf = Pf(:).*H(:);
        end
        
        % store the weighted and unweighted population outputs
        %------------------------------------------------------------------
        layers.unweighted(ins,ij,:) = ( Pf             )     ;% * exp(real(P.L(ins)));
        layers.weighted  (ins,ij,:) = ( Pf * (J(Ji(ij))) );% * exp(real(P.L(ins)));
        layers.iweighted (ins,ij,:) = ( Pf * (J(Ji(ij))) );% * exp(real(P.L(ins)));
        
    end   
end   

%end

clear Pf

% Now compute node proper CSDs from sum of states spectral responses
%--------------------------------------------------------------------------
for inx = 1:ns
    for iny = 1:ns
        if inx ~= iny
            Pf(:,inx,iny) = sum(layers.iweighted(inx,:,:),2) .* conj( ...
                sum(layers.iweighted(iny,:,:),2) );
        else
            if length(Ji) > 1
                Pf(:,inx,iny) = sum((squeeze(layers.iweighted(ins,:,:))),1);%sum(layers.iweighted(inx,:,:),2);
            else
                Pf(:,inx,iny)=squeeze(sum(layers.iweighted(inx,:,:),2));
            end
        end

    end
end

% Addition of system noise & Lead field scaling: L, Gu, Gn, Gs [if requested]
%--------------------------------------------------------------------------
for ins = 1:ns

    Pf0 = Pf(:,ins,ins);

    %if isfield(M,'nosmooth') && M.nosmooth;

    %else
    %   Pf0 = atcm.fun.awinsmooth(Pf0,6);
    %end

    Pf(:,ins,ins) = abs( Pf0(:));


    % Electrode gain: rescale to sum of data spectrum
    %----------------------------------------------------------------------
    SY = sum(spm_vec(M.y));
    Pf(:,ins,ins) = Pf(:,ins,ins) ./ sum(Pf(:,ins,ins) );
    Pf(:,ins,ins) = Pf(:,ins,ins) * SY;

    %SY = max(spm_vec(M.y));
    %Pf(:,ins,ins) = Pf(:,ins,ins) ./ max(Pf(:,ins,ins) );
    %Pf(:,ins,ins) = Pf(:,ins,ins) * SY;

    %Pf(:,ins,ins) = exp(P.L(ins))*Pf(:,ins,ins);


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
else
    y = Pf;
end

s = timeseries;
g = ts;

end

% Helper subfunctions
function R = confounds(w,h)

X0 = gaubasis(length(w),h)';
R  = speye(length(w)) - X0*X0';

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

function g = thresh_fft(x,Ppf,pk,M,ci,ins,w)
pk  = Ppf>(mean(Ppf)+std(Ppf)./x);
Ppf = atcm.fun.makef(w,w(pk)-1,Ppf(pk),ones(length(find(pk)),1)*2);
g = sum( (spm_vec(squeeze(M.y{ci}(:,ins,ins))) - Ppf(:) ).^2 );
end

function [x] = sq(x)
% squeeze function for csds
if size(x,3) > 1, x = squeeze(x); else, x = x(:); end

end

 %Ppf = atcm.fun.match_energy_signals(spm_vec(M.y),Ppf);

            % this is equivalent of: compute the abs fourier transform at FoI
            % %------------------------------------------------------------
            % f      = (1/dt) * (0:(N/2))/N;
            % data   = fd;
            % data   = (data/N);
            % L2     = floor(N/2);
            % data   = data(1:L2+1);
            % S1     = abs(data);
            % w1     = f;

            %Ppf = interp1(w1,full(abs(S1)),w,'linear','extrap') ;
            %Ppf = abs(Ppf);

            % *6.2831853
            %Ppf = atcm.fun.agauss_smooth(Ppf,1);



%for iw = 1:length(w)
            %    wi(iw) = atcm.fun.findthenearest(w(iw),f);
            %end
            %F   = F(:,wi);
            %F   = atcm.fun.HighResMeanFilt(F,1,4);
            %Ppf = (ys*F)./N;  

            %B   = gaubasis(length(w)*2,20);
            %B   = B(:,1:length(w));
            %b   = B'\Ppf(:);
            %Ppf = b'*B; 


            %Ppf = agauss_smooth(Ppf,1);

            %Ppf  = agauss_smooth(abs(Ppf),1);

           % Ppf = atcm.fun.awinsmooth(Ppf,2);
            
            %Ppf = abs(Ppf);
            
            %[Pps]  = atcm.fun.agauss_smooth_mat(Ppf,1.5);         
            %Ppf    = sum(Pps);


            %Ppf    = interp1(f,full(Ppf),w,'linear','extrap') ;%.* (1+w./w(end));
            
            %Ppf = agauss_smooth(Ppf,1);

            %sfun = @(x) 1 - (corr(spm_vec(M.y),agauss_smooth(Ppf,x)').^2);

            %X = fminsearch(sfun,1);

            %Ppf = agauss_smooth(Ppf,1);

            % for the GP created from VtoGauss see
            % https://peterroelants.github.io/posts/gaussian-process-tutorial/

    % residual = spm_vec(M.y) - spm_vec(Pf0);
    % B = gaubasis(length(w),8);
    % b = B'\residual;
    % 
    % resbasis = b.*B;
    % 
    % Pf0 = Pf0(:) + spm_vec(exp(P.d(:)')*resbasis);
    % 

    % Furnish with parameterised Eigenspectrum of innovations {CSD}
    %B   = -[128 64 32 64]'   + 1j*2*pi*[4 12 48 64]' * exp(P.a(2));
    %H   = spm_s2csd(B(1:4),w);
    %Pf0 = Pf0(:).*sum(1000*H,2);

    %Pf0 = atcm.fun.awinsmooth(Pf0,4);
   
    % optimise the smoothness of the vector to match data
    %if ~isfield(M,'dmd')
        %Pf0 = atcm.fun.smooth_optimise(Pf0,M.y{:},1e-1);
    %end


            %N   = length(t);
            %S1  = fd*dt;
            %w1  = ((1:N) - 1)/(N*dt);

            %N         = length(t);
            %S1        = fd*dt;
            %w1        = ((1:N) - 1)/(N*dt);
            %j         = w1 < max(w);
            %S1        = S1(j);
            %w1        = w1(j);

            %j   = w1 < max(w);
            %S1  = S1(j);
            %w1  = w1(j);

            %S1 = atcm.fun.awinsmooth(S1,8);

            %B  = gaubasis(length(w1),round(length(w1)/4));
            %b  = B'\S1';
            %S1 = b'*B; 

            %S1  = agauss_smooth(S1,1);

            % [Pps]  = atcm.fun.agauss_smooth_mat(abs(S1),2);  
            % S1    = sum(Pps);

            %S1 = medfilt1(S1);

            %S1 = envelope(S1,1,'peak');

            %pp = csaps(w1.',S1.');



    % % dfdg
    % Pop = squeeze(layers.unweighted(ins,:,:));
    % J   = exp(P.J);
    % J   = J(find(J));
    % 
    % gfun = @(b) sum( (spm_vec(M.y) - spm_vec(b(:)'*Pop)).^2 );
    % dfdj = jaco(gfun,J,~~J/8,0,1);
    % newJ = J - (dfdj/gfun(J)/8);
    % Pf0 = newJ(:)'*Pop;

                        % % thalamo-cortical delay effects
                        % ct = 8*exp(P.CT); %60;
                        % tc = 3*exp(P.TC); %20; was 3, now 7 Usrey 2000
                        % 
                        % TCDM = [0  0  0  0  0  0  tc tc;
                        %         0  0  0  0  0  0  tc tc;
                        %         0  0  0  0  0  0  tc tc;
                        %         0  0  0  0  0  0  tc tc;
                        %         0  0  0  0  0  0  tc tc;
                        %         0  0  0  0  0  0  tc tc;
                        %         ct ct ct ct ct ct 0  0;
                        %         ct ct ct ct ct ct 0  0];
                        % TCDM = repmat(TCDM,[nk nk]);
                        % v    = v + ((TCDM.*max(dfdx,1e-4))*v);
                        


%             elseif UseSmooth == 2
% 
%                 K = 24;
%                                 
%                 % Compute Dynamic Mode Decomposition on TF data
%                 [Ppf,~,c] = atcm.fun.tfdecomp(pc,dt,w,4,2,@median);
%                 
%                 [q,r]= qr(c);
%                 
%                 sh = hnoisebasis(length(w),exp(P.a));
%                 
%                 Ppf = q*diag(sh)*r;
%                 
%                 %Ppf = smooth(Ppf,22,'loess');
%                 
%                 %Ppf = Ppf(:);%.*Gn(:);
%                 
%                 g = AGenQ(Ppf);
%                 
%                 [Eigenvalues, ev] = atcm.fun.dmd(c, 1, dt);
%                 
%                 %Ppf = Ppf(:).*hamming(length(Ppf));
%                 
%                 %Ppf = Ppf .* hnoisebasis(length(w),exp(P.a));
%                 
%                 %Ppf = ev'*hnoisebasis(K,exp(P.a));
%                 
%                 Ppf = abs( Ppf*ev' );
%                 
%                 
%                 
%                 %[Ppf,~,Pfm] = atcm.fun.AfftSmooth( pc, dw./dt, w, smth) ;
%                 
%                 %SK = AGenQ(std(c,[],2));
%                                 
%                 % Remove ~1/f system noise
%                 %c = abs(c) ./ repmat(Gn(:),[1 size(c,2)]);
%                 
%                 %c = atcm.fun.delembed(pc,K+1)';
%                 
%                 %c = atcm.fun.assa(pc,K+1);
%                 
%                 %[~,~,c] = atcm.fun.AfftSmooth( pc, dw./dt, w, K+1) ;
%                 %c = squeeze(c);
%                 
%                 % Eigenvectors of the DMD of TF matrix are spectra of modes
%                 %warning off;
%                 %[Eigenvalues, ev, ModeAmplitudes, ModeFrequencies, GrowthRates, POD_Mode_Energies] = atcm.fun.dmd(c, K-2, dt);
%                 % warning on;
%                  
% %                  Kern = c*ev';
% %                  
% %                  for ik = 1:K
% %                      Pp(ik,:) = atcm.fun.Afft( Kern(:,ik)', dw./dt, w) ;
% %                      Pp(ik,:) = smooth(Pp(ik,:),8);
% %                  end
% %                  
% %                  Ppf = sum(Pp,1);
% %                  
% %                  nd  = size(P.d,1);
% %                  X   = spm_dctmtx(length(w),nd + 1);
% %                  Mu  = exp(X(:,2:end)*P.d);
% %                  Ppf = Ppf(:).*Mu(:);
%                  
%                 
%                 % Compute modal frequencies
% %                 ModeFrequencies=(angle(Eigenvalues)/pi)*(1/(2*dt));
% %                                 
% %                 % Make values positive (as per a psd estimate)
% %                 ModeAmplitudes  = abs(ModeAmplitudes);
% %                 ModeFrequencies = abs(ModeFrequencies);
% %                 ModeFrequencies = round(ModeFrequencies*100)./100;
% %                 
% %                 %[ModeFrequencies ModeAmplitudes]
% %                 
% %                 % remove duplicates
% %                 [~,I]=unique(ModeFrequencies);
% %                 
% %                 ModeFrequencies = ModeFrequencies(I);
% %                 ModeAmplitudes  = ModeAmplitudes(I);
% %                 
% %                 % Remove stuff outside freqs of interest
% %                 G = find( ModeFrequencies>=w(1) & ModeFrequencies<=w(end) );
% %                                 
% %                 %wt = round(ModeFrequencies)/max(w);
% %                 %ModeAmplitudes = ModeAmplitudes(:).*wt(:);
% %                       
% %                 %Ppf = atcm.fun.makef(w,ModeFrequencies,ModeAmplitudes,ones(length(ModeAmplitudes),1)*1.5);
% %                 
% %                 % Convert to (Gaussian) series\ and average
% %                 for ijf = 1:length(G)
% %                    PF(ijf,:) = atcm.fun.makef(w,ModeFrequencies(G(ijf)),ModeAmplitudes(G(ijf)),1);
% %                 end
% %                  
% %                  % average
% %                  Ppf = spm_vec(max(PF));
% %                  %Ppf = atcm.fun.moving_average(Ppf,2);
%                  %Ppf = AGenQ(Ppf)*Ppf(:);
%                                 
%                 % dct transform
%                 %nd  = size(P.d,1);
%                 %X   = spm_dctmtx(length(w),nd + 1);
%                 %Mu  = exp(X(:,2:end)*P.d);
%                 %Ppf = idct( dct(Ppf(:)).*Mu(:) );
%                                 
%             elseif UseSmooth == 3
%                 
%                 [Ppf,~,c] = atcm.fun.tfdecomp(pc,dt,w,8,   1,@max);
%                 
%                 [g] = AGenQ(std(c'));
%                 
%                 %Ppf = GL*Ppf(:);
%                 
%                 Ppf = g*hnoisebasis(length(w),exp(P.a));
%                 
%                 %Ppf = AGenQ(Ppf)*Ppf;
%                 %Ppf = AGenQ(std(c'))*Ppf;
%                 %Ppf = AGenQ(Ppf)*Ppf;
%                 
%                 nd  = size(P.d,1);
%                 X   = spm_dctmtx(length(w),nd + 1);
%                 Mu  = exp(X(:,2:end)*P.d);
%                 Ppf = idct( dct(Ppf(:)).*Mu(:) );
%                 
%             elseif UseSmooth == 4
%                  
%                 Ppf = atcm.fun.Afft( pc, dw./dt, w)' ;
%                 Ppf = atcm.fun.moving_average(Ppf,2);
%                 
%                 [~,~,c] = atcm.fun.tfdecomp(pc,dt,w,12,2,@mean);
%                 
%                 C = cov(c');
%                                 
%                 % Laplacian smoothing
%                 A  = real(C) .* ~eye(length(C));
%                 N  = size(A,1);
%                 GL = speye(N,N) + (A - spdiags(sum(A,2),0,N,N))/4;
%                 
%                 Ppf = GL*Ppf;
%                 Ppf = GL*Ppf;
%                 Ppf(Ppf<0) = 0;
%                 %Ppf = abs(Ppf);
%                 
%                                 
%                 nd  = size(P.d,1);
%                 X   = spm_dctmtx(length(w),nd + 1);
%                 Mu  = exp(X(:,2:end)*P.d);
%                 Ppf = idct( dct(Ppf(:)).*Mu(:) );

                
                % end -----------------------------------------------------
                
                
                                
%                 %series = atcm.fun.AGenQn(pc,8);
%                 NM     = 2;%12;%rank(series);%12;
%                 
%               %  series = atcm.fun.bandpassfilter(series,1/dt,[w(1) w(end)]);
%                 
%                 % see https://cwrowley.princeton.edu/theses/tu.pdf, pg27
%                 %[Eval, Evec, Amp, Freq, GR, ME] = atcm.fun.dmd(full(cov(series)),NM,dt);
%                 [Eval, Evec, Amp, Freq, GR, ME] = atcm.fun.dmd(full((series)),NM,dt);
%                 
%                 %projection = Evec*series;
%                 %projection = (dfdx*Evec')'*series;
%                 
%                 projection = Evec;
%                 
%                 %NM = 12;%rank(series);
%                 %projection = atcm.fun.assa(sum(series,1),NM)';
%                 
%                 
%                 for i = 1:NM
%                     %Pf(i,:) = atcm.fun.AfftSmooth(projection(i,:),1./dt,w,30);
%                     %Pf(i,:) = atcm.fun.Afft(projection(i,:),1./dt,w);
%                     
%                     if isfield(M,'fooof') && M.fooof
%                        Pf(i,:) = atcm.fun.component_spectrum(w,Pf(i,:),12);
%                     end
%                     
%                     %Pf(i,:) = atcm.fun.awinsmooth(Pf(i,:),4);
%                     %Pf(i,:) = atcm.fun.AGenQn(Pf(i,:),4)*Pf(i,:)';
%                     
%                     %Pf(i,:) = pyulear(projection(i,:),12,w,1/dt);
%                     
%                     [Pf(i,:),~,mat] = atcm.fun.tfdecomp(projection(i,:),dt,w,8,1);
%                     
%                     %Pf(i,:) = max(abs(mat),[],2);
%                 end
%                 
%                 %lmf = M.y{:}\Pf';
%                 %Ppf = lmf*Pf;
% 
%                 %
%                %Ppf = max(Pf,[],1);
%                Ppf = mean(Pf,1);
%               % Ppf = atcm.fun.AGenQn(Ppf,4)*Ppf';
%                 %Ppf = max(Pf);
%              %   Ppf = atcm.fun.awinsmooth(Ppf,8);
%                 Ppf = Ppf(:);
%                 



%                 % Compute modal frequencies
%                 ModeFrequencies=(angle(Eval)/pi)*(1/(2*dt));
%                                 
%                 % Make values positive (as per a psd estimate)
%                 Amp  = abs(Amp);
%                 Freq = abs(Freq);
%                 Freq = round(Freq*100)./100;
%                 
%                 %[ModeFrequencies ModeAmplitudes]
%                 
%                 % remove duplicates
%                 [~,I]=unique(Freq);
%                 
%                 Freq = Freq(I);
%                 Amp  = Amp(I);
%                 
%                 % Remove stuff outside freqs of interest
%                 G = find( Freq>=w(1) & Freq<=w(end) );
%                                 
%                 %wt = round(ModeFrequencies)/max(w);
%                 %ModeAmplitudes = ModeAmplitudes(:).*wt(:);
%                       
%                 %Ppf = atcm.fun.makef(w,ModeFrequencies,ModeAmplitudes,ones(length(ModeAmplitudes),1)*1.5);
%                 
%                 % Convert to (Gaussian) series\ and average
%                 for ijf = 1:length(G)
%                    PF(ijf,:) = atcm.fun.makef(w,Freq(G(ijf)),Amp(G(ijf)),3);
%                 end
%                  
%                 
%                 Ppf = max(PF);
%                 Ppf = Ppf(:);
                                
%                  [Eval, Evec, Amp, Freq, GR, ME] = atcm.fun.dmd(series,NM,dt);
%                  %S = diag(Amp)*Evec;
%                  S = Evec;
%                  NM = atcm.fun.findthenearest(cumsum(abs(Eval))./sum(abs(Eval)),.99);
%                  for i = 1:NM; 
%                      %Pf(i,:) = atcm.fun.AfftSmooth(S(i,:),1./dt,w,20);
%                      
%                      %Pf(i,:) = pyulear(S(i,:),8,w,1/dt);
%                      
%                      Pf(i,:) = atcm.fun.tfdecomp(S(i,:),dt,w,8,2);
%                      
%                      
%                      %Pf(i,:) = atcm.fun.Afft(S(i,:),1./dt,w);
%                      %Pf(i,:) = atcm.fun.AGenQn(Pf(i,:),4)*Pf(i,:)';
%                  end
%                 
%                  %Ppf = max(Pf);
%                  Ppf = mean(Pf,1);
                 
                 %b = atcm.fun.lsqnonneg(atcm.fun.AGenQn(real(Ppf),8),real(Ppf)');
                 %Ppf = atcm.fun.makef(w,w(find(b)),Ppf(find(b)),ones(length(find(b)),1)) ;
                 
                 %Ppf = atcm.fun.awinsmooth(Ppf,4);
                 
                 %Q = atcm.fun.AGenQn(Ppf,8);
                 %Q = Q .* ~eye(length(w));
                 %Q = Q ./ norm(Q);
                 
                 %Ppf = Q*Ppf(:);
                 
                 %Ppf = atcm.fun.AGenQn(Ppf,8)*Ppf(:);
                 
                 
                %signal = Evec*series;
                
                %[Ppf,~,MM] = atcm.fun.AfftSmooth( ifft(Ppf,length(t)), dw./dt, w, 32) ;
                
                %Ppf = atcm.fun.AfftSmooth(pc,dw./dt,w,20)' ;
                
              % X  = spm_dctmtx(nf,9);
              % X  = exp( X(:,2:end)*P.d(1:8) );
                
                %Ppf = pyulear(pc,8,w,1/dt);
                
                %[Ppf,~,MM] = atcm.fun.tfdecomp(Evec*series,dt,w,4,2);
                
                %X = atcm.fun.hnoisebasis(nf,exp(P.a));
                
               % Ppf = Ppf(:);%.*X(:);
                
                %Ppf = pyulear(pc,8,w,1/dt);
                
                % estimate autocov
%                 Ppm  = atcm.fun.estcov(Ppf,length(Ppf)); 
%                 Ppm = Ppm + Ppm';
%                 Ppm = cov(Ppm');
%                 
%                 % comute covariance and smooth with Gauss kernel                
%                 Ppm = Ppm.*atcm.fun.GaussEye(nf)*Ppm';           
%                                 
%                 % eigenvectors are spectral modes s.t. Gaussian
%                 [V,D] = eig(Ppm);
%                 [~,order] = sort(diag(D),'descend');
%                 D = diag(D);
%                 D = D(order);
%                 V = V(:,order); 
%                 
%                 Ppf = spm_vec( V(:,1)'*Ppm );
%                 Ppf = abs(Ppf);
                
                %Ppf = diag(V*V');%*X;
                %Ppf = AGenQ(Ppf)*Ppf;
                

                %                     % no input
%                     if i > 1
%                         % 4-th order Runge-Kutta method.
%                         %--------------------------------------------------
%                         drive(i)=0;
% 
%                         k1 = f(v0          ,drive(i),P,M);
%                         k2 = f(v0+0.5*dt*k1,drive(i)+dt/2,P,M);
%                         k3 = f(v0+0.5*dt*k2,drive(i)+dt/2,P,M);
%                         [k4] = f(v0+    dt*(k3),drive(i),P,M);
%                         
%                         dxdt = (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
%                         v0         = v0 + dxdt;
% 
%                         % State Delays - interpolated
%                         %--------------------------------------------------
%                         d = 100*[.006 .002 .001 .004 .001 .008 .001 .001].*exp(P.ID);
%                         d = repmat(d,[1 nk]);
%                         L = (d);
%                         
%                         for j = 1:length(L)
%                           ti = real(L(j))/dt;
%                           if i > 1 && any(ti)
%                               pt = t(i) - ti;
%                               if pt > 0
%                                 v0(j) = interp1(t(1:i), [y0(j,1:i-1) v0(j)]', pt);
%                               end
%                           end
%                         end
% 
%                         % Full update
%                         %--------------------------------------------------
%                         y0(:,i) =   v0;
% 
%                     else
%                         % 4-th order Runge-Kutta method.
%                         [k1,J] = f(v0          ,drive(i),P,M);
%                         k2 = f(v0+0.5*dt*k1,drive(i),P,M);
%                         k3 = f(v0-0.5*dt*k2,drive(i),P,M);
%                         k4 = f(v0+    dt*k3,drive(i),P,M);
% 
%                         dxdt      = (dt/6)*(k1+2*k2+2*k3+k4);
%                         v0         = v0 + dxdt;
%                         y0(:,i)    = v0 ;
%                     end

                %Ppf = atcm.fun.Afft(pc,1/dt,w);

                %[dPf,in] = atcm.fun.padtimeseries(Ppf);
                %dPf = atcm.fun.HighResMeanFilt(dPf,1,2);
                %Ppf  = dPf(in);

                %[ev,evec] = atcm.fun.dmd(hx,1,dt);
                %Ppf = hx*evec';

                %Ppf = atcm.fun.Afft(pc,1/dt,w);
                %Ppf = atcm.fun.approxlinfitgaussian(real(Ppf),[],[],4)';

                %H   = gradient(gradient(Ppf));
                %K   = 2*exp(P.a(1));
                %Ppf = Ppf - K*H;

                %mu = []; amp = []; wid = [];
                %for i = 1:length(I)
                %    mu  = [mu(:); I{i}.mu(:)];
                %    amp = [amp(:); I{i}.amp(:)];
                %    wid = [wid(:); I{i}.wid(:)];
                %end

                %f = @(pw) atcm.fun.makef(w,mu,amp,wid.*pw);
                %Q = @(x) VtoGauss(spm_vec(x));
                %g = @(pw) norm( (Q(M.y)-Q(f(pw))) * (Q(M.y)-Q(f(pw)))' );

                %inwid = ones(size(wid));
                %[X,F] = fminsearch(g,inwid);

                %Ppf = f(X);

                %H = gradient(gradient(Ppf));
                %Ppf = Ppf - 2*H;

                %Ppf = atcm.fun.AfftSmooth(pc,1/dt,w,20);
              
                %Ppf = pyulear(pc,20,w,1/dt);


                % mu = real(mean(hx,2));
                % st = std(hx,[],2);
                % 
                %  for i = 1:length(mu)
                %    Q(i,:) = atcm.fun.makef(w,i,mu(i),st(i)/2);
                %  end
                %  Q = (Q ./norm(Q))';
                %  [u,s,v] = svd(Q);


                 %V = (Q'+Q)./2;



                

                %V = VtoGauss(real(DCM.xY.y{:}));
                % [q,r] = qr(V');
                % r = abs(r)';
                % I = atcm.fun.findthenearest(cumsum(max(r))./sum(max(r)),.9); I = I(1);
                % r = r(:,1:I);
                % b = find( mean(r) < 0 );
                % r(:,b) = [];
                % 
                % Ppf = mean(r,2);


                % 
                %  [E,D] = eig(Q);
                %  D = diag(D);
                %  [D,I] = sort(D,'descend');
                %  E = E(:,I);
                % 
                %  K = 20;
                %  FDD = E(:,1:K).*D(1:K)';
                %  Ppf = abs(sum((FDD),2));

                 %[u,s,v] = svd(Q);
                 %Q = u(:,1:20)'*Q;
                 %B = Q./norm(Q);
                
                %Ppf = ProcessSignal(w(:),real(Ppf),0,3,3,0,0,0,0,0,0);;
                
                
                
                %[Ppf,I] = gausscoeffit(Ppf,20);                
                %B = gaubasis(length(w),8);
                %B=sqrt(B);
                %b = B'\real(Ppf);
                %Ppf = abs(b'*B);
               

                % Convert to (~Gaussian) features matrix 
                %Ppf = atcm.fun.approxlinfitgaussian(real(Ppf),[],[],2);
