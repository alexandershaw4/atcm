function [y,w,s,g,t,pst,layers,noise,firing,QD] = integrate7(P,M,U,varargin)
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
% replacement for spm_csd.mtf. 
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
%    g   :    J'* y                      (the weighted signal / LFP timeseries)
%    k[i]:    L * ( J[i] * fft(y[i]) )   (over populations, i) 
%  c/psd :   (k * Gu) + Gn(auto) + Gs(cross)       
%
% Fourier transforms are computed using AfftSmooth.m, which computes cross
% spectra using FFT and the complex-conjugate method. Smoothing is done by
% chunking the series into overlapping windows, computing fft and
% averaging.
%
% Noise model - 3 components (Gu, Gn, Gs)
% -------------------------------------------------------------------------
%   Gu(:,i) = exp(P.a(1,i))*(w.^0);                
%   Gn      = exp(P.b(1) - 2)*w.^(-exp(P.b(2)))'; 
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



% FoI & initial states
%--------------------------------------------------------------------------
w    = M.Hz;
x    = M.x;

% See if dt | fs is specified & generate sampletimes & inputs
%--------------------------------------------------------------------------
dt    = 1/1200;                   % hard-wired 400 hz
Fs    = 1/dt;                     % sampling frequency
tn    = 4;%2;%1.6;                % sample window length, in seconds
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


% Check motion: f(x) & integration: IS(f(x))
%--------------------------------------------------------------------------
f   = spm_funcheck(M.f); 
x   = atcm.fun.solvefixedpoint7(P,M);
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
    [y{c},w,s{c},g{c},t{c},layers{c},noise{c},firing{c},QD{c}] = dodxdt(pst,f,v,Q,M,dt,w,drive);

end



end

function [y,w,s,g,t,layers,noise,firing,QD] = dodxdt(t,f,v,P,M,dt,w,drive)
% Numerical integration, signal weighting and FFT of timeseries with
% spline interpolation and smoothing


% Choose either the DCM kernels scheme, or a Euler update scheme
%--------------------------------------------------------------------------
warning off;
IntMethod  = ' ';
WithDelays = 2;   % 45 = 4th order Runge-Kutta method with delays

% Prerequisits for integration with delays
if WithDelays == 2 || WithDelays == 5
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

switch IntMethod

    case 'kernels'
        
        % Just use the kernels approach?
        %-----------------------------------
        [K0,K1,K2,H1] = spm_kernels(M,P,length(t),dt);
        y             = H1(:,1:length(v))';
        S = [];

    otherwise
                
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
                    v = v + Q*f(v,drive(i),P,M);
                end        
                y(:,i) = v;
                       
                
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
            VR     = -40;
            Vx     = exp(P.S)*32;
            V      = spm_unvec(v,M.x);
            S(:,:,i) = spm_Ncdf_jdw(V(:,:,1),VR,Vx);     % mean firing rate  
            
        end
end

warning on;

% Reshape to model state space outputs
%--------------------------------------------------------------------------
[ns,npp,nk] = size(M.x);
y(isnan(y)) = 0;
y = reshape(y,[ns npp nk size(y,2)]);
timeseries  = y;
firing      = S;
nf          = length(w);

% Options
%--------------------------------------------------------------------------
%DoSmoothing    = 1; % smooth the fourier series with a 10% Loess window
DoHilbert      = 0; % take the absolute (magnitude) of the hilbert envelope
Bandpassfilter = 0; % band pass filter to [w(1)-1) : w]
DoDCT          = 0; % discrete cosine transform series before fft
IncDCS         = 0; % include a cosine set in the noise model


% Compute channel noise before computing spectral response
%--------------------------------------------------------------------------
% Neuronal innovations: a multiplier on the model signal
for i = 1:size(P.a,2)
    %Gu(:,i) =  exp(P.a(1,i))*(w.^(-exp(P.a(2,i))));
    Gu(:,i) = exp(P.a(1,i))*(w.^0);
end

% Spectrum of channel noise (non-specific): added to spectrum
Gn = 0*P.b(1)*(w.^0)';

% Spectrum of channel noise (specific): added to spectrum
for i = 1:size(P.c,2)
    %Gs(:,i) = exp(P.c(1,i) - 2)+w.^(-exp(P.c(2,1)));
    Gs(:,i) = exp(P.c(1,i) )+w.^(-exp(P.c(2,1)));
    %Gs(:,i) = exp(P.c(1,i) - 2)*(w.^0);
    Gs = Gs*0;
end

% % Neuronal innovations: a discrete cosine set
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
end

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
burn = findthenearest(500,t); 

% % Generate weighted (principal cells) signal - g & fft(g)i
% %--------------------------------------------------------------------------
J      = ( exp(P.J(:)) ) ;
Ji     = find(J);
ts     = zeros(ns,length(t));

for ins = 1:ns
    
    % J-weighted sum of this channel / region
    %----------------------------------------------------------------------
    xseries   = full( squeeze(y(ins,:,:,:)) );
    xseries   = reshape(xseries, [npp*nk,length(t)] );
    ts0       = (J'*xseries)';
        
    % transform to (smoothed) cosine set if requested
    %----------------------------------------------------------------------
    if DoDCT
        ts0 = dct(ts0);
    end
        
    % apply electrode gain & store this channel / node 
    %----------------------------------------------------------------------
    ts(ins,:) = ts(ins,:) + ts0'  * exp(P.L(ins));
    
end

% Compute transforms of individual cell / layer signals
%--------------------------------------------------------------------------
J           = ( exp(P.J(:)) ) ;
Ji          = find(J);
J           = reshape(J,[npp nk]);
[ijx,ijy,z] = find(J);                      % indices of active (~0) states
pfYl        = zeros(ns,length(w));

for ins = 1:ns
    for ij = 1:length(Ji)
        y0 = squeeze(y(ins,ijx(ij),ijy(ij),:))'; % region n, state ij
        
        if Bandpassfilter
            for i = 1:size(y0,1)
                chunk   = y0(i,:);
                [chy,ii]= atcm.fun.padtimeseries(chunk);                   % pad fun
                %chy     = atcm.fun.bandpassfilter(chy, 1/dt,[w(1)-3 100]); % filt fun
                chy     = lowpass(chy,4,1/dt);
                y0(i,:) = chy(ii);
            end
        end    
        
        % Spectra of this state, with innovations
        %------------------------------------------------------------------
        [Pf,Hz]  = atcm.fun.AfftSmooth(y0(burn:end),1/dt,w); %*
        
        for i = 1:length(Hz)
            Pf(i,:,:) = sq(Pf(i,:,:))*diag(Gu(i,:))*sq(Pf(i,:,:))';
        end
        
        % Smoothing
        %------------------------------------------------------------------
        Pf       = smooth(Hz,Pf,0.2,'loess');
        Pf       = full(Pf);
        J        = full(J);
        
        % Actively contributing states, weighted
        %------------------------------------------------------------------
        layers.weighted(ins,ij,:)   = abs( Pf * J(ijx(ij),ijy(ij)) ) ;%/ length(w) ;        
        layers.weighted(ins,ij,:)   = layers.weighted(ins,ij,:) * exp(P.L(ins));
        pfYl(ins,:)                 = pfYl(ins,:) + squeeze(layers.weighted(ins,ij,:))';
        
        % Actively contributing states, weighted (with imaginary, for CSD)
        %------------------------------------------------------------------
        layers.iweighted(ins,ij,:)   =   ( Pf * J(ijx(ij),ijy(ij)) ) ;%/ length(w) ;        
        layers.iweighted(ins,ij,:)   = layers.iweighted(ins,ij,:) * exp(P.L(ins));
        pfYl(ins,:)                  = pfYl(ins,:) + squeeze(layers.weighted(ins,ij,:))';
        
        % Actively contributing states, unweighted
        %------------------------------------------------------------------
        layers.unweighted(ins,ij,:) = ( Pf ) * exp(P.L(ins));
        
    end
        
end

% Now compute CSD from sum of weighted cells per region
%--------------------------------------------------------------------------
for inx = 1:ns
    for iny = 1:ns
        if ~DoHilbert
            Pf(:,inx,iny) = sum(layers.iweighted(inx,:,:),2) .* conj( ...
                            sum(layers.iweighted(iny,:,:),2) );
        else
            Pf(:,inx,iny) = max(hilbert(layers.iweighted(inx,:,:)),[],2) .* conj( ...
                            max(hilbert(layers.iweighted(iny,:,:)),[],2) );

        end
    end
end

% Take the absolute (magnitude) of the cross spectra
%--------------------------------------------------------------------------
Pf = abs(Pf)/length(w);

% Incorporate noise components for auto (Gs) and cross (Gn) spectra
%--------------------------------------------------------------------------
for i = 1:ns
    Pf(:,i,i) = Pf(:,i,i) + Gs(:,i);
    for j = 1:ns
        Pf(:,i,j) = Pf(:,i,j) + Gn;
    end
end

% returns for this trial - {g}
%--------------------------------------------------------------------------
y = real(Pf);
w = Hz;
s = timeseries;
g = ts;
t = drive;

end

function [x] = sq(x)
% squeeze function for csds
if size(x,3) > 1, x = squeeze(x); else, x = x(:); end

end
