function [y,w,s,g,t,pst,layers,noise,firing,QD,Spike] = integrate_linearise(P,M,U,varargin)
%
% AS


% w, initial states, dt | fs is specified & generate sampletimes & inputs
%--------------------------------------------------------------------------
w     = M.Hz;                     % FoI (w)
x     = M.x;                      % model (hidden) states
Kx    = x;                        % pre-fp x, for kernel scheme
dt    = 1/1200;                   % hard-wired 1200 hz
Fs    = 1/dt;                     % sampling frequency
tn    = 2;                        % sample window length, in seconds
pst   = 1000*((0:dt:tn-dt)');     % peristim times we'll sample

% unpack simulation pst options, if specified
%--------------------------------------------------------------------------
if isfield(M,'sim')
    %tn  = M.sim.pst(end);
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
        
    case 4
        
        % TWO oscillatory inputs...
        %------------------------------------------------------------------
        mu1   = .001*exp(P.R(1));                      % mean amplitude
        mu2   = .001*exp(P.R(2));
        mf1   = 50*exp(P.R(3));                  % frequency
        mf2   = 10*exp(P.R(4));
        
        drive(:,2) = .2*(mu1 * sin(2*pi*mf1*(pst/1000)) ) ;  % (sin) oscillation over time
        drive(:,1) = .2*(mu2 * sin(2*pi*mf2*(pst/1000)) );
        
        
end


% expansion (fixed) point: trial & parameter effects are deviations from here
%--------------------------------------------------------------------------
f    = spm_funcheck(M.f); 
if solvefp; x    = atcm.fun.solvefixedpoint(P,M);
else ;      x    = x;
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


% Initialise states [v] and period [t] for integration loop
%--------------------------------------------------------------------------
M.pst = t;

% Linearise & first order Taylor expansion
%------------------------------------------------------

% Here we aim to linearise f, such that:
% dx = f(x) --> dx = Ax
%
% then we use matrix A to predict the timeseries of the system

% initial point
x0      = spm_vec(M.x);

% find a fixed point
xbar    = atcm.fun.solvefixedpoint(P,M);
M.x     = spm_unvec(xbar,M.x);

% compute Jacobian, evaluated at xbar
dfdx    = spm_diff(M.f,M.x,M.u,P,M,1);

% eigenvectors and values
[T,D]   = eig(full(real(dfdx)));
iT      = pinv(T);

% integrate: x(t) = T*exp(D*t)*iT*x0
for i = 1:length(t)
    y(:,i) = T*exp(D*(t(i)./1000))*iT*x0;
end
S=[];


% Reshape to model state space outputs
%--------------------------------------------------------------------------
[ns,npp,nk] = size(M.x);
y(isnan(y)) = 0;
y(isinf(y)) = 0;
y           = reshape(y,[ns npp nk size(y,2)]);
timeseries  = y;
firing      = S;
nf          = length(w);
spike       = [];

if isfield(M,'IncDCS')
    IncDCS = M.IncDCS;
end



% Neuronal innovations: a multiplier on the model signal
%--------------------------------------------------------------------------
for i = 1:size(P.a,2)
    %Gu(:,i) =  exp(P.a(1,i))*(w.^(-exp(P.a(2,i))));
    Gu(:,i) = exp(P.a(1,i))*(w.^0);                  % P.a = constant
end

% Spectrum of channel noise (non-specific): added to spectrum
%--------------------------------------------------------------------------
Gn = P.b(1)*(w.^0)';                                 % P.b = constant

% Spectrum of channel noise (specific): added to spectrum
%--------------------------------------------------------------------------
for i = 1:size(P.c,2)
    Gs(:,i) = exp(P.c(1,i) )+w.^(-exp(P.c(2,1)));    % P.c = expone
end

if HamNoise
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
    Gu = repmat( (exp(P.a(1))*Mu) , [1 ns]);
end

% Return noise components for this trial
%--------------------------------------------------------------------------
noise.Gu = Gu;
noise.Gn = Gn;
noise.Gs = Gs;


% Implement the observation model:  y = [ L * fft(J' * y) ] + noise
%--------------------------------------------------------------------------    
% Optional burn in - i.e. take transform of n:end ms instead of 0:end...
burn = findthenearest(300,t); 

% Generate weighted (principal cells) signal - g & fft(g)i
%--------------------------------------------------------------------------
J      = ( exp(P.J(:)) ) ;
Ji     = find(J);
ts     = zeros(ns,length(t));

% This is the main loop over regions to calculate the spectrum
%--------------------------------------------------------------------------
for ins = 1:ns
        
    % extract time series of all states from this region
    %----------------------------------------------------------------------
    yx = reshape( squeeze(y(ins,:,:,:)), [npp*nk,length(t)]); 
    yx = yx(1:size(y,2),:);        % limit to membrane potentials 
    Eigenvectors = yx;

    % compute spectra for populations & reduce (fit) using Gaussian mixtures
    %----------------------------------------------------------------------
    for ij = 1:length(Ji)
        
        y0 = Eigenvectors(ij,:) ;
        Hz = w;
        
        % just a smoothed fft of the (contributing) states
        [Pf,Hz]  = atcm.fun.AfftSmooth(Eigenvectors(Ji(ij),burn:end),1/dt,w,20); 
        Pf = (abs(Pf));
        %Pf = smooth( Pf , 8 ,'lowess' );
        
        % reduce each to a unimodela distribution
        w1 = w./w(end); Pf = Pf.*w1(:);
        [curve] = fit((w).',Pf,'gauss3');
        
        c(1,:) = curve.a1*exp(-((w-curve.b1)/curve.c1).^2)';
        c(2,:) = curve.a2*exp(-((w-curve.b2)/curve.c2).^2)';
        c(3,:) = curve.a3*exp(-((w-curve.b3)/curve.c3).^2)';
        
        Pf = sum(c,1)';
        Pf = smooth( Pf , 8 ,'lowess' );
        
        Pf0(ins,ij,:) = Pf;
    end      
end

clear Pf

% Do Dynamic Mode Decomposition on the frequency space response
%--------------------------------------------------------------------------
layers.iweighted = [];
layers.weighted = [];
for ins = 1:ns
    
    Pfa   = squeeze( (Pf0(ins,:,:)));
        
    X1 = Pfa(:,1:end-1);
    X2 = Pfa(:,2:end);
    
    % Dynamic mode decomposition
    [U,S,V] = spm_svd(X1);
    
    r  = length(P.pca);
    U  = U(:,1:r);
    S  = S(1:r,1:r);
    V  = V(:,1:r);
    
    Atilde    = full(U')*X2*full(V)*pinv(full(S));
    [W,eigvs] = eig(Atilde);
    Phi       = X2*V*pinv(full(S))*W;
    Q         = real(Phi'*Pfa)';
    
    layers.phi(ins,:,:) = Phi';
    
    for im = 1:length(P.pca)
        WC(im,:) = exp(P.pca(im))*(Q(:,im))';
        
        if DoHamming
            WC(im,:) = WC(im,:).*H(:)';
        end
        if isfield(P,'psmooth')
            warning off
            WC(im,:) = smooth( WC(im,:) , (P.psmooth(im)) ,'lowess' );
            warning on
        end
        
        %Multiply in the semi-stochastic neuronal fluctuations
        WC(im,:) = WC(im,:).*Gu(:,ins)';
        
    end
    
    % store
    layers.iweighted(ins,:,:) = WC * exp(P.L(ins));
    layers.weighted = layers.iweighted;
    
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
for i = 1:ns
    %Pf(:,i,i) = Pf(:,i,i) + Gs(:,i);
    for j = 1:ns
        Pf(:,i,j) = Pf(:,i,j) + Gn;                % TURN ON!!!!!!!
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

% % kill final portion of frequency window - if data wand bandpass filtered
% %--------------------------------------------------------------------------
% if KillTail
%     % bandpass filteing the real data forces the edges to ~0, but the model
%     % can't necessarily do that, so force it here 
%     kf   = 70; % tail off from 70
%     H    = ones(nf,1); H(end) = 0;
%     Ikf  = atcm.fun.findthenearest(w,kf); 
%     Wexp = 1./(1:length(Ikf:nf));
%     Wexp = smooth(Wexp,.3);
%     for i = 1:ns
%         for j = 1:ns
%             Pf(Ikf:end,i,j) = Pf(Ikf:end,i,j).*Wexp;
%         end
%     end
% end

% returns for this trial - {g}
%--------------------------------------------------------------------------
y = Pf;
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