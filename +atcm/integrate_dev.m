function [y,w,s,g,t,pst,layers,other] = integrate_dev(P,M,U,varargin)
% Development version of integrate_1

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
        mf    = 20*exp(P.R(2));                      % frequency
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
        mu1   = 0.5*exp(P.R(1));                % mean amplitude
        mu2   = 1.0*exp(P.R(2));
        mf1   = 50*exp(P.R(3));                  % frequency
        mf2   = 10*exp(P.R(4));
        
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
end

if ~isfield(M,'timefreq')
    M.timefreq = 0;
end

% expansion (fixed) point: trial & parameter effects are deviations from here
%--------------------------------------------------------------------------
f    = spm_funcheck(M.f); 

% solve for a fixed point, or not
if solvefp; x    = atcm.fun.solvefixedpoint(P,M);
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
p     = abs(eigs(dfdx,1,'SR',OPT));
N     = ceil(max(1,dt*p*2));
n     = spm_length(M.x);
dfdx  = dfdx - speye(n,n)*exp(-16);
Q     = (spm_expm(dt*D*dfdx/N) - speye(n,n))/dfdx;
%Q    = sum( ((D.^N).*J)*(Q^N)/factorial(N) ):
QD    = Q;        

% Set up custom delay operator - parameterised state x state delays
% conferred by population delays (c.f. conduction delay)
%--------------------------------------------------------------------------
if isfield(P,'ID')
    
    if npp == 8 
        del = exp(P.ID).*[2 1/4 1/2 8 1/2 4 2 2]/2.4;
        del = exp(P.ID).*[4 1/4 1 1/2 1/2 4 2 1]/2.4;
        del = exp(P.ID).*[1 1 1/2 1 1/2 1 1 1];
                
        if isfield(P,'IDb')
            deli = exp(P.IDb).*[1 1/5 1/2 1 1/2 1 1/2 8]/2.4;
            deli = exp(P.IDb).*[4 1/5 1 1/2 1/2 4 2 144]/2.4;
        end
        
        del = repmat(del,[1 nk]);
        del = 1./del;
        
        if isfield(P,'IDb')
            deli = repmat(deli,[1 nk]);
            deli = 1./deli;
        else
            deli=del;
        end
        
        if ns > 1
            if ~isfield(P,'delay')
                del = (spm_vec(repmat(del,[ns 1])))';
            else
            end
        end
        condel=del;
        
    elseif npp == 2
        del = exp(P.ID).*[1 1]./2.4;
        del = repmat(del,[1 nk]);
        del = 1./del;
        deli=del;
    end
    
else
    del    = 1;
    condel = 1;
end

% convert parameterised delay vector to state-by-state matrix
%--------------------------------------------------------------------------
if isfield(P,'ID')
    del = sqrt(del)'*sqrt(deli);
else
    del = 1;
end

% Can also pass a full np-by-np delay matrix (that would be a lot of params)
if isfield(P,'ID') && all(size(P.ID)==8)
    del = [4 1/4 1 8 1/2 4 2 20]/2.4;
    del = exp(P.ID).*( sqrt(del)'*sqrt(del) );
    del = repmat(del,[nk nk]);
    del = 1./del;    
end

% Recompute operator with del incorporated - i.e. approximate a full DDE
% using a delay operator (Q) following Q ~ mexp(dt*D*del*dfdx./N - I)
%--------------------------------------------------------------------------
if isfield(P,'ID')
    deldfdx = del.*dfdx;
    deldfdx(isnan(deldfdx))=0;
    deldfdx(isinf(deldfdx))=0;
    p       = abs(eigs(deldfdx,1,'SR',OPT));
    N       = ceil(max(1,dt*p*1));
    n       = spm_length(M.x(:));
    dfdx    = dfdx - speye(n,n)*exp(-16);
    %Q       = (spm_expm(dt*D*(del.*dfdx)/N) - speye(n,n))/(dfdx);
    Q       = (spm_expm(dt*D*(deldfdx)/N) - speye(n,n))/(dfdx);
end
 
% Use a different delay operator (within population {diagonal}) with the
% Newton-Cotes integration algorithm
%--------------------------------------------------------------------------
if WithDelays == 44
    if isfield(P,'ID')
        if npp == 8
            del = exp(P.ID).*[.01 1.2 1 1 1 1 .08 .08]; % (tau = 1./x)
        elseif npp == 2
            del = exp(P.ID).*[1.2 1];
        end
    else
        del = [1 1];
    end
    
    del = repmat(del,[1 nk]);
    del = 1./del;
    %Q = spm_expm(dt*diag(del)*dfdx/(N*n)); % matrix exponential diagonal delay operator
    Q = spm_expm(dt*diag(del));
end

condel = 1;

% firing rate & count when fired (membrane potential passes threshold)
%--------------------------------------------------------------------------
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

switch IntMethod

    case 'kernels'
        
        % Just use the kernels approach?
        %------------------------------------------------------------------
        [K0,K1,K2,H1] = spm_kernels(M,P,length(t),dt);
        y = H1(:,1:end-1)';
        S = [];
        
    case 'ode45'
        % matlab build in ode solver
        ode = @(t,v,P,M,f) spm_vec( Q*f(spm_unvec(v,M.x),drive,P,M) );
        [~,y]   = ode113(ode,t/1000,spm_vec(v),drive,P,M,f);
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
        
        % initial point
        x0      = spm_vec(M.x);
        
        % find a fixed point
        if solvefp; xbar    = spm_vec(atcm.fun.solvefixedpoint(P,M));
        else;       xbar    = spm_vec(M.x);
        end
        M.x     = spm_unvec(xbar,M.x);
        
        % compute (states) Jacobian, evaluated at xbar
        dfdx    = spm_diff(M.f,M.x,M.u,P,M,1);
            
        % compute du/dt (inputs)
        dfdu    = spm_diff(M.f,M.x,drive(1),P,M,2);        
        
        % eigenvectors and values
        [T,D]   = eig(full(real(del.*dfdx)));
        iT      = pinv(T);
        d       = diag(D);
                
        % integrate: x(t) = T*exp(D*t)*iT*x0 
        for i = 1:length(t)
            % We still want to drive this (linear) system - i.e.
            % x(t) = T*exp(D+(dfdu*input)*t)*iT*x0 
            
            x_inp  = dt*dfdu*drive(i);
            y(:,i) = T*diag(exp((d+x_inp)*(t(i)./1000)))*iT*xbar;
            %y(:,i) = T*diag(exp(d*(t(i)./1000)))*iT*xbar;
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

            elseif WithDelays == 2 
                % Karl's Euler-like-with-a-matrix exponential delay operator
                % this just essentially an RK method
                for j = 1:N
                    %v  = v + del'.*(Q*f(spm_unvec(v,M.x),drive(i,:),P,M));
                    %v0 = v0 + del'.*(Q*f(spm_unvec(v0,M.x),0.0001,P,M));                     
                    
                    %v  = v + ((del.*Q)*f(spm_unvec(v,M.x),drive(i,:),P,M));
                    %v0 = v0 + (del.*Q*f(spm_unvec(v0,M.x),0.0001,P,M)); 

                    v  = v + (Q*f(spm_unvec(v,M.x),drive(i,:),P,M));
                    v0 = v0 + (Q*f(spm_unvec(v0,M.x),0.0001,P,M)); 
                end
                y(:,i) = v; %- spm_vec(M.x);  
                y0(:,i) = v0;
                
                
            elseif WithDelays == 44
                % Newton-Cotes quadrature, including Q
                for j = 1:N
                    k1 = f(v      ,drive(i,:),P,M);
                    k2 = f(v+dt*k1,drive(i,:),P,M);
                    k3 = f(v+dt*k2,drive(i,:),P,M);
                    k4 = f(v+dt*k3,drive(i,:),P,M);
                    k5 = f(v+dt*k4,drive(i,:),P,M);

                    dxdt = ((2*dt)./45)*(7*k1 + 32*k2 + 12*k3 + 32*k4 + 7*k5);
                    v    = v + Q*dxdt;
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

% Break phase-locking if modelling an induced response
RandPhase = 1;
if RandPhase
    for i = 1:size(y,1)
        tmp = fft(y(i,:));
        y(i,:) = ifft(real(tmp)+sqrt(-1)*fliplr(imag(tmp)));
    end
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
try
    yy = spm_unvec(y0,y);
catch
    yy = y;
end

series.States_without = yy;
series.States_with_inp = y;

% Compute the cross spectral responses from the integrated states timeseries 
%==========================================================================

% EVOKED spectrum: fft(f(input) - f(no_input))
%------------------------------------------------------
%[y,s,g,noise,layers] = spectral_response(P,M,y-yy,w,npp,nk,ns,t,nf,timeseries,dt,dfdx,ci,1);

% INDUCED spectrum: fft(f(input)) - fft(f(no_input))
%------------------------------------------------------
% Randomise phase if modelling an induced response
% RandPhase = 1;
% if RandPhase
%     for j = 1:10;
%         for i = 1:size(y,1)
%             tmp = fft(y(i,:));
%             r   = rand(size(tmp));
%             r   = exp( 2*pi*1i*r );
%             yph(i,:) = ifft(tmp.*r);
%         end
%         [y1{j},s,g,noise,layers1] = spectral_response(P,M,yph,w,npp,nk,ns,t,nf,timeseries,dt,dfdx,ci,1);
%     end
% end
% 
% y1 = mean(cat(2,y1{:}),2);

% system spectral response with specified input
[y1,s,g,noise,layers1] = spectral_response(P,M,y,w,npp,nk,ns,t,nf,timeseries,dt,dfdx,ci,1);
%[y0,s,g,noise,layers0] = spectral_response(P,M,yy,w,npp,nk,ns,t,nf,timeseries,dt,dfdx,ci,1);

%y = spm_unvec( (spm_vec(y1) - spm_vec(y0))./spm_vec(y0), y1);
%y = spm_unvec( (spm_vec(y1) - spm_vec(y0)), y1);

y = y1;
%y = spm_unvec( real(spm_vec(y)), y);
%y(y<0) = -y(y<0);

series.with_inp = y;

%y = full(exp(P.Ly)*y);
t = drive;

layers = layers1;
%layers = spm_unvec( (spm_vec(layers1)-spm_vec(layers0)), layers1);

end

function [y,s,g,noise,layers]=spectral_response(P,M,y,w,npp,nk,ns,t,nf,timeseries,dt,dfdx,ci,type)
% Main spectral response function with lots of options.

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
        
        y0 = Eigenvectors(ij,:) ;
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
                test = this;
                pc = test;
                nn = size(test,1);
            end
            
            % Remove super slow drifts
            pc = detrend(pc);
            
            clear Ppf  Pfm Ppf1 Ppf2 Ppf3

            % Select sliding-window smoothed fft or whole chunk fft
            %--------------------------------------------------------------
            UseSmooth = 1;
            if UseSmooth
                
                % User specified FFT smoothing (num windows)
                %----------------------------------------------------------
                if isfield(M,'smth') && ~isempty(M.smth)
                    smth = M.smth;
                else
                    smth = 30;
                end

                % Smoothed Fourier transform (atcm.fun.AfftSmooth)
                %----------------------------------------------------------
                for i = 1:nn; 
                    [Ppf(i,:),~,Pfm(i,:,:)] = atcm.fun.AfftSmooth( pc, dw./dt, w, smth) ;
                end

                % Retain principal eigenmode over (Gaussian) time-freq windows
                Ppf = squeeze(Pfm)';
                [u,s,v] = svd(Ppf');n = 1;
                Ppf = spm_vec(u(:,n)*s(n,n)*mean(v(:,n)));
    
            else % (else use non smooth)
            
                % Non-smoothed fft version:
                %--------------------------------------------------------------
                for i = 1:nn; Ppf = atcm.fun.Afft( pc(i,:), dw./dt, w) ;end

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
            
            % Distribution mixture model - estimate distributions & system noise
            %--------------------------------------------------------------
            [~,pk] = atcm.fun.maxpoints(Pf,50); % peak n points
            wint   = 1:length(w);
            
            % Dist fits are on residuals - ie data minus already explained.
            % This is to avoid all populations converging on same portion of spectrum 
            %--------------------------------------------------------------
            if ij == 1; ResY = squeeze(M.y{ci}(:,ins,ins));
            else;       ResY = ResY - Pf0;
            end
            
            % This function compares fitting with Gaussian, Cauchy, Laplace 
            % and Gamma dists and returns best fit
            %-------------------------------------------------------------- 
            Pf0  = atcm.fun.findbestdist(wint,wint(pk)-1,Pf(pk),2*ones(length(find(pk)),1),ResY);
            Pf   = Pf0;
            
        end
    
        % Make it non-sparse and vectorised
        warning on;
        Pf = (Pf)';
        Pf = full(Pf)';
        J  = full(exp(P.J));
        
        if DoHamming
            H = .5+hamming(nf,'periodic');
            [val,ind] = max(H);
            H(1:ind)  = val;
            Pf = Pf(:).*H(:);
        end

        % store the weighted and unweighted population outputs
        %------------------------------------------------------------------
        layers.unweighted(ins,ij,:) = ( Pf             )      * exp(P.L(ins));
        layers.weighted  (ins,ij,:) = ( Pf * abs(J(Ji(ij))) ) * exp(P.L(ins));
        layers.iweighted (ins,ij,:) = ( Pf * abs(J(Ji(ij))) ) * exp(P.L(ins));
        
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
Pf = (Pf)/length(w);


% Addition of system noise & Leadfield scaling: L, Gu, Gn, Gs
%--------------------------------------------------------------------------
for ins = 1:ns
    
    % Electrode gain
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
