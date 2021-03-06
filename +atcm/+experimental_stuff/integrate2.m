function [y,w,s,g,t,pst,layers,noise,firing,QD,Spike] = integrate2(P,M,U,varargin)
% Numerical integration and spectral response of a neural mass model.
% This is the NETWORK version, which handles networks of connected models,
% and with different trial types.
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
% and neuronal innovations (discrete cosine set of order = number of
% opopulations) are computed and returned.
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
% Observation model:
% -------------------------------------------------------------------------
%    g    :    J'* y + noise              (the weighted signal / LFP timeseries)
%    Pf[i]:    L * ( J[i] * fft(y[i]) )   (over populations, i) 
%  c/psd  :    (Pf * Gu) + Gn(auto) + Gs(cross)       
%
% Fourier transforms are computed using AfftSmooth.m, which computes cross
% spectra using FFT and the complex-conjugate method. Smoothing is done by
% chunking the series into overlapping windows, computing fft and
% averaging.
%
% Noise model - 3 components (Gu, Gn, Gs)
% -------------------------------------------------------------------------
%   Gu(:,i) = exp(P.a(1,i))*(w.^0);                
%   Gn      = 0;
%   Gs(:,i) = exp(P.c(1,i) - 2)*w.^(-exp(P.c(2,1)));
%
% Neuronal Innovations
% -------------------------------------------------------------------------
%   dct = sqrt(2/N)*cos(pi*(2*n+1)*(k-1)/(2*N))
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
    
    Q       = (spm_expm(dt*D*dfdx/N) - speye(n,n))*spm_inv(dfdx);
    QD      = Q;
    
    %Q = D*dt;

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
% if solvefp; M.x    = atcm.fun.solvefixedpoint(P,M);
% else ;      M.x    = M.x;%*0;
% end
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
        
        % Also compute layer-specific spectra
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
        
    case 'aMTF'
        
        
%         % Alex's MTF
%         j    = adfdx(f,P,M,1); j(isnan(j))=0;
%         dfdu = adfdx(f,P,M,2); dfdu(isnan(dfdu))=0;
%         dgdx = spm_diff(M.g,M.x,M.u,P,M,1);
%         j    = Q*j;
%         dfdu = Q*dfdu;
%         
%         %[v0,s0] = eig(j);
%         %s0 = diag(s0);
%         
%         [u0,s0] = eig(j);
%         s0 = diag(s0);
%         %y  = u0*sin( exp(s0)*t' );
%         %1./(1j*2*pi*w - s(k));
%         s0 = 1j*imag(s0) + real(s0) - exp(real(s0));
                
%         t0 = t/t(end);
%         y  = pinv(u0)*sin( 2*pi*exp(s0*t0') );
%         
%         y  = pinv(u0)*sin( 2*pi*exp(s0*t0') );
%         
%         %y = sin(j*v*t');
%         
%         S = [];

    otherwise
        
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
                    
                    % Ozaki 1992 numerical method:
                    % A bridge between nonlinear time-series models and
                    % nonlinear stochastic dynamical systems: A local 
                    % linearization approach:
                    % "dx(t) = (expm(dfdx*t) - I)*inv(dfdx)*f"
                    %[fx,dfdx] = f(v,drive(i),P,M);
                    %v = v + spm_dx(D*dfdx,D*fx,dt);
                                        
                end   
                % Expansion point - i.e. deviation around fixed point
                y(:,i) = v - spm_vec(M.x);
                
%             elseif WithDelays == 3
%                 % stepwise updates of the delay matrix: recompute jacobian
%                 
%                 M0   = M;
%                 M0.x = spm_unvec(v,M.x);
%                 j    = adfdx(f,P,M0,2); j(isnan(j))=0;
%                 OPT.tol = 1e-6*norm((j),'inf');
%                 if OPT.tol == 0;OPT.tol = 1;end
%                 
%                 p    = abs(eigs(j,1,'SR',OPT));
%                 N    = ceil(max(1,dt*p*2));
%                 N    = min([N 4]);
%                 n    = spm_length(M.x);
%                 Q    = (spm_expm(dt*D*j/N) - speye(n,n))*spm_inv(j);
%                 QD   = Q;
%                 
%                 for j = 1:N
%                     %v = v + Q*f(v,drive(i),P,M,Curfire);
%                     v = v + Q*f(v,drive(i),P,M);           % CHANGE ME BACK
%                 end   
%                 % Expansion point - i.e. deviation around fixed point
%                 y(:,i) = v - spm_vec(M.x);
                
                
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
%                 VR       = -40;
%                 Vx       = exp(P.S)*32;
%                 V        = spm_unvec(v,M.x);
%                 Curfire  = spm_Ncdf_jdw(V(:,:,1),VR,Vx);     % mean firing rate  
%                 S(:,:,i) = Curfire;
                S = [];
                 %Dfire   = (1./[1 1/2 1/4 1 1/2 1 8 8]);
                %Dfire   = (1./[1 1   1/4 1 1/2 1 1 1]);   % CHANGE ME BACK!!
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

%fprintf('finished\n');

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
    yw          = dy;
    WithDelays  = 20; % flag to invoke fft(yw)
    y = yy;
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


% % phase-orthogonalised signal noise
%--------------------------------------------------------------------------
% y0     = reshape(y,[ns*npp*nk,length(t)]);
% %[ph,PhsCor] = atcm.fun.computephase({y});
% ph = y0;
% ph     = ph{1};
% ph     = symm_orthog( reshape(ph,[ns*npp*nk,length(t)]) );
% ph     = reshape(ph,[ns*npp*nk,length(t)]);
% D      = ( (-D)/dt )*10;
% D0     = D;
% yn     = y0*0;
% for it = 1:length(t)
%     for i = 1:size(v,1)
%         for j = 1:size(v,2)
%             try
%                 ph_diff  = ph(i,it-D0(i,j)) - ph(j,it);
%                 yn(i,it) = yn(i,it) + ph_diff;
%             end
%         end
%     end
% end



% Spectral Response Options
%--------------------------------------------------------------------------
%DoSmoothing    = 1; % smooth the fourier series with a 10% Loess window
DoHilbert      = 0; % take the absolute (magnitude) of the hilbert envelope
Bandpassfilter = 0; % band pass filter to [w(1)-1) : w]
DoDCT          = 0; % discrete cosine transform series before fft
IncDCS         = 1; % include a cosine set in the noise model                % off
DoHamming      = 1; % Cosine Hann window the model spectrum      
HamNoise       = 0; % Hann noise components - if data was BPF, exponential 
                    % delay based noise model wont fit without hanning edges
KillTail   = 0;
DoPCA      = 0;
DoSpecResp = 1; % 1=fft(integrated signal), 2=fft(volterra kernels)

% Compute channel noise before computing spectral response
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

%[Gu,Gs,Gn] = spm_csd_mtf_gu(P,M.Hz);

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
    %if size(Mu,2) == 1, Mu = Mu*ones(1,ns); end        % UNCOMMENT ME!!!!
    %Gu = Gu.*Mu;
    Gu = exp(P.a(1))*Mu;
end


% Return noise components for this trial
%--------------------------------------------------------------------------
noise.Gu = Gu;
noise.Gn = Gn;
noise.Gs = Gs;


% Implement observation model:  y = [ L * fft(J' * y) ] + noise
%--------------------------------------------------------------------------
% Computes the (smoothed) FFT of the weighted (J) signal (g = J*y).
% Adds parameterised channel noise and neuronal innovations in the form of
% a exponential decay functions and a discrete cosine set of order = number
% of populations, respectively.
%
% Use a lengthy burn-in to ensure the system (& optimsation) are in
% steady-state (generating oscillations).
%
    
% Optional burn in - i.e. take transform of n:end ms instead of 0:end...
burn = findthenearest(300,t); 
%burn = findthenearest(0  ,t); 

% Generate weighted (principal cells) signal - g & fft(g)i
%--------------------------------------------------------------------------
J      = ( exp(P.J(:)) ) ;
Ji     = find(J);
ts     = zeros(ns,length(t));

if WithDelays == 20 && ~strcmp(IntMethod,'kernels')
    % this method has already calculated the weighted signal
    ts = yw;
        
    % append some noise
    s2n = mean(ts(:))*1/8;
    ts  = ts + s2n*randn(size(ts));
    %ts  = ts + 1/8*shuffle(ts);

else
    % otherwise compute it now
    for ins = 1:ns

        % J-weighted sum of this channel / region
        %----------------------------------------------------------------------
        xseries   = full( squeeze(y(ins,:,:,:)) );
        xseries   = reshape(xseries, [npp*nk,length(t)] );
        ts0       = (~~J'*xseries)';           % use only J>0 states
        %ts0 = (sparse(1:8,1,1,size(xseries,1),1)'*xseries)'; % use all mV states

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
end

if DoPCA == 1
    yy  = squeeze(y(1,:,1,:));
    yy         = PEig(yy',1);
    y(1,:,1,:) = yy';
end

% for inx = 1:ns
%     % remove a few modes from states time series
%     x = squeeze( y(ins,:,:,:) );
%     x = reshape( x , [npp*nk, length(t)] );
% 
%     [Eigenvalues, Eigenvectors, ModeAmplitudes, ModeFrequencies, ...
%         GrowthRates, POD_Mode_Energies]=dmd_rom(x, 4, dt);
%     y(ins,:,:,:) = reshape(Eigenvectors,[npp,nk,length(t)]);
% end

% for inx = 1:ns
%     % remove a few modes from states time series
%     x = squeeze( y(ins,:,:,:) );
%     x = reshape( x , [npp*nk, length(t)] );
%     
%     [u0,s0,v0] = spm_svd(x);
%     p1 = u0(:,1)*s0(1,1)*v0(:,1)';
%     %p2 = u0(:,2)*s0(2,2)*v0(:,2)';
%     x = x - p1 ;%- p2;
%     y(ins,:,:,:) = reshape(x,[npp,nk,length(t)]);;
% end
%timeseries=y;

% remove prinicipal eigenmode from mV series
% x = squeeze( y(1,:,1,:) );
% [u0,s0,v0] = spm_svd(x);
% p1 = u0(:,1)*s0(1,1)*v0(:,1)';
% p2 = u0(:,2)*s0(2,2)*v0(:,2)';
% x = x - p1 - p2;
% y(1,:,1,:) = x;
% timeseries=y;


if DoSpecResp == 1
    
    % Compute transforms of individual cell / layer signals
    %----------------------------------------------------------------------
    J           = ( exp(P.J(:)) ) ;             % contirbuting states
    Ji          = find(J);                      % CS indices
    J           = reshape(J,[npp nk]);          % 'state' space matrix
    [ijx,ijy,z] = find(J);                      % indices of active (~0) states
    pfYl        = zeros(ns,length(w));

    % normal - with lots of components...1 else fft only=0
    SpecNoise    = 0;
    
    if SpecNoise == 1
        for ins = 1:ns
            for ij = 1:length(Ji)
                y0 = squeeze(y(ins,ijx(ij),ijy(ij),:))'; % region n, state i,j
                %y0 = y0-mean(y0);
                %y0 = detrend(y0);

                if Bandpassfilter                        % ill advised
                    for i = 1:size(y0,1)
                        chunk   = y0(i,:);
                        [chy,ii]= atcm.fun.padtimeseries(chunk);                    % pad fun
                        %chy     = atcm.fun.bandpassfilter(chy, 1/dt,[w(1)-3 100]); % filt fun
                        chy     = lowpass(chy,4,1/dt);
                        y0(i,:) = chy(ii);
                    end
                end    

                % Spectra of this state, with innovations
                %--------------------------------------------------------------
                [Pf,Hz]  = atcm.fun.Afft(y0(burn:end),1/dt,w); %*
                Pf = Pf';

                for i = 1:length(Hz)
                    Pf(i,:,:) = sq(Pf(i,:,:))*diag(Gu(i,:))*sq(Pf(i,:,:))';
                end

                % Smoothing
                %--------------------------------------------------------------
                Pf = full(Pf);
                J  = full(J);

                if ~DoPCA
                    %Pf = full(atcm.fun.HighResMeanFilt(Pf,1,4));
                    Pf = smooth(Hz,Pf,0.20,'loess');   
                end

                if DoHamming
                    H  = kaiser(nf,2.5);
                    H(1:round(nf/2)) = 1;
                    Pf = Pf.*H;
                end

                % Actively contributing states, weighted
                %--------------------------------------------------------------
                layers.weighted(ins,ij,:)   = abs( Pf * J(ijx(ij),ijy(ij)) ) ;%/ length(w) ;        
                layers.weighted(ins,ij,:)   = layers.weighted(ins,ij,:) * exp(P.L(ins));
                pfYl(ins,:)                 = pfYl(ins,:) + squeeze(layers.weighted(ins,ij,:))';

                % Actively contributing states, weighted (with imaginary, for CSD)
                %--------------------------------------------------------------
                layers.iweighted(ins,ij,:)   =   ( Pf * J(ijx(ij),ijy(ij)) ) ;%/ length(w) ;        
                layers.iweighted(ins,ij,:)   = layers.iweighted(ins,ij,:) * exp(P.L(ins));
                pfYl(ins,:)                  = pfYl(ins,:) + squeeze(layers.weighted(ins,ij,:))';

                % Actively contributing states, unweighted
                %--------------------------------------------------------------
                layers.unweighted(ins,ij,:) = ( Pf ) * exp(P.L(ins));

            end

        end
    
    elseif SpecNoise == 0 
        % do a straight fft, no noise!
        for ins = 1:ns
            for ij = 1:length(Ji)
                y0 = squeeze(y(ins,ijx(ij),ijy(ij),:))'; % region n, state i,j
                       
               %[Pf,Hz]  = atcm.fun.Afft(y0(burn:end),1/dt,w);
               %Pf = Pf';
               
               %Pf = sgolayfilt(Pf, 5, length(w))'; 
               %spf = SmoothFAS(Pf',.1,17);
               %Pf = [Pf(1); spf];
               
               % uncomment
               warning off;
               
               try
                   [Pf,Hz]  = atcm.fun.AfftSmooth(y0(burn:end),1/dt,w,60); 
               catch
                    Hz = w;
                    Pf = w'*0;
               end
                
               warning on;
               
               
               
                %[Pf,Hz] = AfftSmoothWin(y0(burn:end),1/dt,w,20);
                %Hz=Hz';
                
                %[Pf,Hz]  = atcm.fun.AfftSmooth(y0 ,1/dt,w,50); 
                
%                 mar = spm_mar(y0(burn:end)',8); 
%                 mar = spm_mar_spectra(mar,w,1/dt);
%                 Pf  = mar.P;
%                 Hz  = w;
%                 [u s v] = spm_svd(Pf,1);
%                 Pr      = u(:,1)*s(1,1)*mean(v(:,1));
%                 Pf      = full(spm_unvec(Pr,Pf));
                
                %[Pf,Hz] = pwelch(y0(burn:end),200,[],w,1/dt);
                %[Pf,Hz] = pcov(y0(burn:end),3,w,1/dt);
                
                
                Pf = Pf';
                %Hz = Hz';
                Pf = full(Pf)';
                Pf = Pf .* Hz';                                            % off
                
                %Pf = smooth(Hz,Pf,0.3,'loess');
                
                %Pf = smooth(Hz,Pf,0.2,'loess');
                
                %s0 =        smooth(       Pf ,0.3,'loess');
                %s1 = flipud(smooth(flipud(Pf),0.3,'loess'));
                
                % Moving average where the smoothing kernel gets bigger
                % as you slide along the vector, from 0 : 30%
%                 ik = round( linspace( 1 , round(.3*length(Pf)) , 6 ) );
%                 for is = 1:length(Pf)
%                     kern = ik( sum(~(is < ik)) );
%                     if is < length(Pf)-ik(end)
%                         Pfz(is) = mean(Pf(is:is+kern));
%                     else
%                         Pfz(is) = mean(Pf(is:end));
%                     end
%                 end
%                Pf = Pfz';
                
                %Pf = smooth(Hz,Pf,0.4,'loess');
                %Pf = smooth(Hz,Pf,0.52,'loess');
                
                % new
                for i = 1:length(Hz)
                    Pf(i,:,:) = sq(Pf(i,:,:))*diag(Gu(i,:))*sq(Pf(i,:,:))';
                end
                
                
                %Pf = smooth(Hz,Pf,13,'sgolay'); 
                %Pf = NewMeanFilt3D(Pf,2,2);
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
                
                % Actively contributing states, weighted
                %--------------------------------------------------------------
                layers.weighted(ins,ij,:)   = abs( Pf * J(ijx(ij),ijy(ij)) ) ;%/ length(w) ;
                layers.weighted(ins,ij,:)   = layers.weighted(ins,ij,:) * exp(P.L(ins));
                pfYl(ins,:)                 = pfYl(ins,:) + squeeze(layers.weighted(ins,ij,:))';

                % Actively contributing states, weighted (with imaginary, for CSD)
                %--------------------------------------------------------------
                layers.iweighted(ins,ij,:)   =   ( Pf * J(ijx(ij),ijy(ij)) ) ;%/ length(w) ;
                layers.iweighted(ins,ij,:)   = layers.iweighted(ins,ij,:) * exp(P.L(ins));
                %pfYl(ins,:)                  = pfYl(ins,:) + squeeze(layers.weighted(ins,ij,:))';

                % Actively contributing states, unweighted
                %--------------------------------------------------------------
                layers.unweighted(ins,ij,:) = ( Pf ) * exp(P.L(ins));
            end
        end
        
        
    end
    

    if WithDelays ~= 20
        
%         for inx = 1:ns
%             % remove principal mode
%             x = squeeze( layers.iweighted(ins,:,:) );
%             [u0,s0,v0] = spm_svd(x);
%             p1 = u0(:,1)*s0(1,1)*v0(:,1)';
%             x = x - p1;
%             layers.iweighted(ins,:,:) = x;
%         end
        
        % Now compute node CSDs from sum of weighted cells per region
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
    else
        
        % if the weighted signal was precomputed, just use an FFT of this
        %------------------------------------------------------------------
        % Smoothed fft of pre-weighted
        Pf = atcm.fun.AfftSmooth(yw(:,burn:end),1/dt,w,50); %*
        Pf = Pf.*w';
        %Pf = atcm.fun.Afft(yw(:,burn:end),1/dt,w)'; %*
        
%         yi        = yw(:,burn:end);
%         Nw        = 1./dt;
%         
%         S1        = fft(yi)'*dt;
%         S1        = S1((1:Nw/2) + 1);
%         w1        = ((1:Nw) - 1)/(Nw*dt);
%     
%         j         = w1 < max(w);
%         Pf        = S1(j,1,1);
%         w1        = w1(j);
% 
%         % crop un-requested lower freqs
%         Pf(1:(w(1)-1)) = [];
%     
%         %Pf = atcm.fun.Afft(yw,1/dt,Hz)';
        
        
        % incorporate Gu
        for i = 1:length(Hz)
            Pf(i,:,:) = sq(Pf(i,:,:))*diag(Gu(i,:))*sq(Pf(i,:,:))';
        end
        
        % smooth
        %Pf = full(atcm.fun.HighResMeanFilt(Pf,1,4));
        %Pf = smooth(Hz,Pf,0.20,'loess');  
    end

    % Take the absolute (magnitude) of the cross spectra
    %----------------------------------------------------------------------
    Pf = abs(Pf)/length(w);

%    if WithDelays ~= 20
        % Incorporate noise components for auto (Gs) and cross (Gn) spectra
        %----------------------------------------------------------------------
        for i = 1:ns
           %Pf(:,i,i) = Pf(:,i,i) + Gs(:,i);
           for j = 1:ns
               Pf(:,i,j) = Pf(:,i,j) + Gn;
               
%                DoNoiseGamma = 1;
%                if DoNoiseGamma
%                    Gam = atcm.fun.makef(w,exp(P.G(1))*60,exp(P.G(2))*.05,2);
%                    Pf(:,i,j) = Pf(:,i,j) + Gam(:);
%                end

               
           end
        end
        
%     elseif WithDelays == 20
%         % if == 20, use spm_csd_mtf style noise
%         [Gu,Gs,Gn] = spm_csd_mtf_gu(P,M.Hz);    
%         G     = zeros(nf,ns,ns);
%         for i = 1:nf
%             G(i,:,:) = sq(Pf(i,:,:))*diag(Gu(i,:))*sq(Pf(i,:,:))';
%         end
%         for i = 1:ns
%             G(:,i,i) = G(:,i,i) + Gs(:,i);
%             for j = 1:ns
%                 G(:,i,j) = G(:,i,j) + Gn;
%             end
%         end
%         
%     end
    
elseif DoSpecResp == 0 % if not DoSpecResp, using kernels
    layers = [];
    Pf(:,1,1) = S1;
    
elseif DoSpecResp == 2
    
    % compute volterra kernels
    %----------------------------------------------------------------------
    % reset expansion point
    M.x    = Kx;
    % compute kernels
    Pf     = spm_csd_mtf(P,M,U);
    Pf     = Pf{1};
    % and update noise for consistency with csd_mtf
    [Gu,Gs,Gn] = spm_csd_mtf_gu(P,M.Hz);
    noise.Gu = Gu;
    noise.Gn = Gn;
    noise.Gs = Gs;

    % layer kernel ffts
    for i = 1:length(Ji)
        P0   = P;
        P0.J = P.J*0 - 1000;
        P0.J(Ji(i)) = P.J(Ji(i));
        layers.weighted(1,i,:) = spm_dcm_mtf(P0,M);
    end
    
end

if strcmp(IntMethod,'kernels')
    layers.weighted(1,:,:) = Sj;
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
y = real(Pf);
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
%                     y(ip,i)    = y(ip,i)    - squeeze( y(ip   ,i-(Df(ip)-1) ) );
%                     %y(ip+8,i)  = y(ip+8,i)  + squeeze( y(ip+8 ,i-(Df(ip)-1) ) );
%                     %y(ip+16,i) = y(ip+16,i) + squeeze( y(ip+16,i-(Df(ip)-1) ) );
%                     %y(ip+24,i) = y(ip+24,i) + squeeze( y(ip+24,i-(Df(ip)-1) ) );
%                     %y(ip+32,i) = y(ip+32,i) + squeeze( y(ip+32,i-(Df(ip)-1) ) );
%                 end
%             end
            

