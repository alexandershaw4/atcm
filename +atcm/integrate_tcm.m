function [y,w,s,g,drive,pst,layers,other] = integrate_tcm(P,M,U,varargin)
% Numerical integration and spectral response of a neural mass model.
%
% This version uses a 4th order Runge-Kutta numerical integration with
% matrix delay operator (using a local linearisation / linearisation
% operator at each step), and the Laplace transform for the spectral
% response.
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
% The RK4 DDE implementation is as follows:
% 
%  k1 = Q*f(v          ,drive(i)     ,P,M);
%  k2 = Q*f(v+0.5*dt*k1,drive(i)+dt/2,P,M);
%  k3 = Q*f(v+0.5*dt*k2,drive(i)+dt/2,P,M);
%  k4 = Q*f(v+    dt*k3,drive(i)     ,P,M);
% 
%  dxdt = (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
%  
%  dxdt = dxdt / dt; % step in seconds;
% 
%  compute a linear update operator for this step
%  b    = pinv(full(J)'.*v).*dxdt;
%  Q    =  (J.*b); % dxdt = Q*x;
% 
%  add delays and recompute step @ dt
%  Q    = Q + (D.*~~J);
%  dxdt = (dt*Q*v);
%  v    = v + dxdt;
%
% At each freq of interest, w(k), Laplace is given by;
%
%  Fs(k) = sum(y .* exp(-s(k) * t'/1000) * dt);
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
        
        % and an ERP like input for exogenous input to thalamic
        fre = 40 * exp(P.R(2));
        drive = exp(P.R(1)) * sin(2*pi*fre*pst./1000);
                  
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

% expansion (fixed) point: trial & parameter effects are deviations from here
%--------------------------------------------------------------------------
f    = spm_funcheck(M.f); 

% solve for a fixed point, or not
if solvefp;
    x = spm_unvec(atcm.fun.alexfixed(P,M,1e-6),M.x);
end

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
        dodxdt(pst,f,v,Q,M,dt,w,drive,Kx,U,method,solvefp,c);

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
                            dodxdt(t,f,v,P,M,dt,w,drive,Kx,U,method,solvefp,ci)

% Numerical integration, signal weighting and FFT of timeseries

warning off;
[ns,npp,nk] = size(M.x);

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
[v,J] = (M.f(M.x,0,P,M));
v0    = v;

% Frequency steps: dw
dw = 1./(w(2)-w(1));

% setup a baseline jacobian using finite differecnes
fun = @(v) f(spm_unvec(v,M.x),drive(1),P,M);


% Do numerical integration for a discrete epoc
%------------------------------------------------------------------
for i   = 1:length(t)

    % Matrix delays and rate (integration) constants
    %--------------------------------------------------
    if i == 1
        [k1,J,D] = f(v    ,0*drive(:,i),P,M);

        D   = real(D);
        ddt = dt;
        J   = full(J);
    end

    % integrate w 4-th order Runge-Kutta method.
    %--------------------------------------------------
    k1   = f(v             ,0*drive(i),P,M);
    k2   = f(v+0.5.*ddt.*k1,0*drive(i)+0.5*dt,P,M);
    k3   = f(v+0.5.*ddt.*k2,0*drive(i)+0.5*dt,P,M);
    k4   = f(v+     ddt.*k3,0*drive(i)+dt,P,M);

    dxdt = (ddt/6).*(k1 + 2*k2 + 2*k3 + k4);

    % work in seconds;
    dxdt = dxdt / dt;

    % compute a linear update operator for this step
    b    = pinv(full(J)'.*v).*dxdt;
    Q    =  (J.*b); % dxdt = Q*x;

    % add delays and recompute step @ dt
    Q    = Q + D;
    dxdt = (dt*Q*v);
    v    = v + dxdt;

    % Full update
    %--------------------------------------------------
    y(:,i) =   (v );


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

warning on;

try
      S = m;
catch S = [];
end

% Reshape to model state space outputs
%--------------------------------------------------------------------------
[ns,npp,nk] = size(M.x);
y(isnan(y)) = 0;
y(isinf(y)) = 0;

% recompute dfdx now we've propoerly burned in 
[fx, dfdx]  = f(spm_unvec(y(:,end),M.x),0,P,M);
yorig       = y;
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
g = [];

% Compute the cross spectral responses from the integrated states timeseries 
%==========================================================================
[yb,s,~,layers] = spectral_response(P,M,y,w,npp,nk,ns,t,nf,timeseries,dt,J,ci,1,drive);
y = yb;

end


function [y,s,g,layers]=spectral_response(P,M,y,w,npp,nk,ns,t,nf,timeseries,dt,dfdx,ci,type,drive)
% Main spectral response function with lots of options.

% Spectral Response Options
%--------------------------------------------------------------------------
DoHamming      = 0; %(1)% Cosine Hann window the model spectrum      

if isfield(M,'DoHamming')
    DoHamming = M.DoHamming;
end

% Implement the observation model:  y = [ L * fft(J' * y) ] + noise
%--------------------------------------------------------------------------
    
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

% adjustment for when frequency intervals ~= 1 Hz; i.e. dt = dw/dt
%--------------------------------------------------------------------------
dw = 1./(w(2)-w(1));

% This is the main loop over regions to calculate the spectrum
%==========================================================================
for ins = 1:ns
        
    % extract time series of all states from this region
    %----------------------------------------------------------------------
    if ndims(y) == 4
        yx = reshape( squeeze(y(ins,:,:,burn:burno)), [npp*nk,length(t(burn:burno))]);
    else
        yx = y(:,burn:burno);
    end
    
    % Spectral responses of (observed) states
    %----------------------------------------------------------------------
    for ij = 1:length(Ji)
        Hz = w;
            
        try    
            clear Ppf

            ys  = yx(Ji(ij),:);
            ys  = atcm.fun.bandpassfilter(ys,1/dt,[w(1) w(end)+1/2]);
            %ys  = real(ys);
            %ys  = ys - mean(ys);

            % % filtering
            ti = t(burn:burno);
            ti = ti - ti(1);


            % Laplace transform (of timeseries)
            %---------------------------------------------------
            sigma = 0 ;              % Real part (decay factor)
            f     = 2*pi*w;
            s     = sigma + 1i * f;  % s values (complex)
            Fs    = zeros(size(s));

             for k = 1:length(s)                    
                 Fs(k) = sum(ys .* exp(-s(k) * ti'/1000) * dt);
             end

            %H  = gradient(gradient(Fs));
            %a  = exp(P.a(1))*1;
            %Fs = Fs + a*H; 

            Ppf = abs(Fs);

            %Ppf = atcm.fun.awinsmooth(Ppf,8);
            %Ppf = atcm.fun.agauss_smooth(Ppf,1);

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
            y  = Pf;
            s  = timeseries;
            g  = ts;
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
                Pf(:,inx,iny) = squeeze(sum(layers.iweighted(inx,:,:),2));
            end
        end

    end
end

% Addition of system noise & Lead field scaling: L, Gu, Gn, Gs [if requested]
%--------------------------------------------------------------------------
for ins = 1:ns

    Pf0 = Pf(:,ins,ins);

    %Pf0 = atcm.fun.agauss_smooth(Pf0,1,[],'laplace');
    
    Pf(:,ins,ins) = ( Pf0(:));

    % Electrode gain: rescale to sum of data spectrum
    %----------------------------------------------------------------------
    Pf(:,ins,ins) = exp(P.L(ins))*Pf(:,ins,ins);


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
% if isfield(P,'f')
%     nw = length(w);
%     f     = (1:nw)';
%     f     = exp(real(P.f(1)) + real((P.f(2))*f)/nw);     
%     for i = 1:ns
%         for j = 1:ns
%             Pf(:,i,j) = Pf(:,i,j).*f;
%         end
%     end
% end

% returns for this trial - {g}
%--------------------------------------------------------------------------
if ns == 1 
    y = (Pf); 
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


function [x] = sq(x)
% squeeze function for csds
if size(x,3) > 1, x = squeeze(x); else, x = x(:); end

end


            %ul = [w(1)*exp(P.f(1)) (w(end)-1)*exp(P.f(2))];
            %ys = conv(ys,abandpassfilter(1/dt,ul,ul,length(ti)),'same')*dt;

            %[Ppf] = pwelch(ys,[],[],w);
            %Ppf = atcm.fun.agauss_smooth(Ppf,1);


                        %dxdt = (ddt/6).*(Q.*D_dt)*v;

                        %b    = pinv(full(J)'.*v).*dxdt;
                        %Q    = (D_dt.*J).*b; % dxdt = Q*x; Q is linear operator
                        %dxdt = (Q*v);

                        % Di = D/dt;
                        % d = b*0;
                        % for ix = 1:56
                        %     for iy = 1:56
                        %         if any(Di(ix,iy)) && i > Di(ix,iy)
                        %             d(ix,iy) = interp1(t(i-1:i), full([J(ix,iy) b(ix,iy)]), t(i-Di(ix,iy)))';
                        %         end
                        %     end
                        % end


   % SY = sum(spm_vec(M.y));
   % Pf(:,ins,ins) = Pf(:,ins,ins) ./ sum(Pf(:,ins,ins) );
   % Pf(:,ins,ins) = Pf(:,ins,ins) * SY;

    %SY = max(spm_vec(M.y));
    %Pf(:,ins,ins) = Pf(:,ins,ins) ./ max(Pf(:,ins,ins) );
    %Pf(:,ins,ins) = Pf(:,ins,ins) * SY;

    % noise;



    %H   = gradient(gradient(Pf0));
    %Pf0 = Pf0 - (exp(P.d(1))*3)*H;
    %smoothfun = @(d) sum((spm_vec(M.y) - spm_vec(atcm.fun.agauss_smooth(Pf0 - (d(1)*3)*H,(1.6)*d(2)))).^2);

    %[X,F] = fminsearch(smoothfun,[1 1]);

    %Pf0 = spm_vec(atcm.fun.agauss_smooth(Pf0 - (X(1)*3)*H,(1.6)*X(2)));

    %[parts,moments]=iterate_gauss(Pf0,2);

    %PfL = squeeze(layers.iweighted);

    %b = PfL'\M.y{:}';%atcm.fun.lsqnonneg(PfL',M.y{:}');
    %b = atcm.fun.lsqnonneg(PfL',M.y{:}');;

    %Pf0 = b'*PfL;

    %[eval,evec] = atcm.fun.dmd(PfL',1,1);

    %Pf0 = evec*PfL;

    %[u,s,v] = svd(PfL);

    %Pf1 = u(1:3,:)*s*v';

    %Pf0 = Pf1(1,:) * exp(P.d(1)) + Pf1(2,:) * exp(P.d(2)) + Pf1(3,:) * exp(P.d(3));

    %b   = PfL'\M.y{:}';
    %5Pf0 = b'*PfL;

    %H   = gradient(gradient(Pf0));
    %Pf0 = Pf0 - (exp(P.d(1))*3)*H;

    %Pf0 = atcm.fun.agauss_smooth(Pf0,1);


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
