function [y,w,s,g,t,pst,layers,noise,firing,QD,Spike] = integrate_kern(P,M,U,varargin)
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
% (1. ) Straight Euler scheme (because delays are in equations, f):
%
%    y(i+1) = (i) + dt*f(x,...P)
%
% (2. ) Euler-like scheme with delays computed from Jacobian [Q=D*J]:
%
%    y(i+1) = y(i) + Q*dt*f(y(i),..,P)
%
% (3. ) Or 4th Order Runge-Kutta with or without delays:
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
dt    = 1/1200;                   % hard-wired 400 hz
Fs    = 1/dt;                     % sampling frequency
tn    = 2;                        % sample window length, in seconds
pst   = 1000*((0:dt:tn-dt)');     % peristim times we'll sample

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

% now override main output with spm_csd_mtf!
[y,w] = spm_csd_mtf(P,M,U);

% copmute layers responses
J = exp(P.J);
Ji = find(J);
layers = [];

for ij = 1:length(Ji)
    QP   = P;
    QP.J = QP.J*0-1000;
    QP.J(Ji(ij)) = J(Ji(ij));
    [Pfj,w] = spm_csd_mtf(QP,M,U);
    layers{1}.weighted(1,ij,:,:) = spm_vec(Pfj);
end

% Check motion: f(x) & integration: IS(f(x))
%--------------------------------------------------------------------------
f   = spm_funcheck(M.f); 
x   = atcm.fun.solvefixedpoint(P,M);
M.x = x;
v   = spm_vec(x);
    

% Integration and spectral response for this trial
%--------------------------------------------------------------------------
for  c = 1:size(U.X,1)
    
    % condition-specific parameters
    %----------------------------------------------------------------------
    Q   = spm_gen_Q(P,U.X(c,:));
        
    % integration, spectral response, firing etc. (subfunction)
    %----------------------------------------------------------------------
    [~,~,s{c},g{c},t{c},~,noise{c},firing{c},QD{c},Spike{c}] = dodxdt(pst,f,v,Q,M,dt,w,drive);
    
end


end

function [y,w,s,g,t,layers,noise,firing,QD,spike] = dodxdt(t,f,v,P,M,dt,w,drive)
% Numerical integration, signal weighting and FFT of timeseries with
% spline interpolation and smoothing

% Integration options:
% WithDelays == 0  : Euler without delays
%            == 2  : Euler with delays
%            == 5  : Runge-Kutta 45 w/ DCM delays
%            == 45 : Runge-Kutta 45 with delays
%

% Choose either the DCM kernels scheme, or a Euler update scheme
%--------------------------------------------------------------------------
warning off;
IntMethod  = ' ';
WithDelays = 20;   % 45 = 4th order Runge-Kutta method with delays
[ns,npp,nk] = size(M.x);

% Prerequisits for integration with delays
if WithDelays == 2 || WithDelays == 5 || WithDelays == 20
    [fx, dfdx,D] = f(M.x,0,P,M);
    OPT.tol = 1e-6*norm((dfdx),'inf');
    p       = abs(eigs(dfdx,1,'SR',OPT));
    N       = ceil(max(1,dt*p*2));
    n       = spm_length(M.x);
    Q       = (spm_expm(dt*D*dfdx/N) - speye(n,n))*spm_inv(dfdx);
    QD      = Q;
else
    [fx, dfdx,Q] = f(M.x,0,P,M);
    QD           = Q;
end

% initial firing rate
Curfire = zeros(size(M.x,2),ns)';
firings = cell(ns,1);
DoSpecResp = 1;

% rediscover expansion point
%M.x = atcm.fun.solvefixedpoint(P,M);

switch IntMethod

    case 'kernels'
        
        % Just use the kernels approach?
        %------------------------------------------------------------------
        [K0,K1,K2,H1] = spm_kernels(M,P,length(t),dt);
        y             = H1(:,1:length(v))';
        S             = [];
        
    otherwise
        
        % Do an actual numerical integration, if not using kernel approach
        %------------------------------------------------------------------
        for i   = 1:length(t)
                        
            if ~WithDelays
                % Use a Euler integration scheme: without Delays
                % y(i+1) = y(i) + dt*f(x,P)
                dxdt   = Q*f(v,drive(i),P,M);
                v      = v + dt*dxdt;                
                y(:,i) = v;
                                               
            elseif WithDelays == 2
                % Karl's Euler-like-with-a-Jacobian-Delay scheme
                % dx = (expm(dt*J) - I)*inv(J)*f(x,u)
                for j = 1:N
                    %v = v + Q*f(v,drive(i),P,M,Curfire);
                    v = v + Q*f(v,drive(i),P,M);
                end   
                y(:,i) = v;
                
            elseif WithDelays == 20
                for j = 1:N
                    %v = v + Q*f(v,drive(i),P,M,Curfire);
                    v = v + Q*f(v,drive(i),P,M);
                end
                
                % Weighted at each integration step!
                y (:,i) = v;
                yw(:,i) = spm_gx_erp(spm_vec(v),drive(i)',P,M);
                                     
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
                    v         = v + dxdt;
                    y(:,i)    = Q*v;       
                    
            end
                
            % firing function at dxdt - assumes conductance model
            %--------------------------------------------------------------
            VR       = -40;
            Vx       = exp(P.S)*32;
            V        = spm_unvec(v,M.x);
            Curfire  = spm_Ncdf_jdw(V(:,:,1),VR,Vx);     % mean firing rate  
            S(:,:,i) = Curfire;
            
            if size(M.pE.H,1) == 8
                % thalamo cortical model
                Dfire   = (1./[1 1/2 1/4 1 1/2 1 8 8]);
                Dfire   = Dfire .* exp(P.TV);
            elseif size(M.pE.H,1) == 6
                % cortex only model
                Dfire   = (1./[1 1/2 1/4 1 1/2 1]);
                Dfire   = Dfire .* exp(P.TV);
            end
            
            Curfire = Curfire .* repmat(Dfire,[ns 1]);
            
            for ins = 1:ns
                fired          = find(squeeze(V(ins,:,1)) >= VR*Curfire); 
                firings{ins}   = [firings{ins}; [i+0*fired',fired'] ];
            end
        end
end

warning on;

% Reshape to model state space outputs
%--------------------------------------------------------------------------
[ns,npp,nk] = size(M.x);
y(isnan(y)) = 0;
y           = reshape(y,[ns npp nk size(y,2)]);
timeseries  = y;
firing      = S;
nf          = length(w);
spike       = firings;

% returns for this trial - {g}
%--------------------------------------------------------------------------
y = [];
%w = Hz;
s = timeseries;
g = yw;
t = drive;
layers = [];
noise = [];
end

function [x] = sq(x)
% squeeze function for csds
if size(x,3) > 1, x = squeeze(x); else, x = x(:); end

end


