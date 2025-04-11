function [Y,w,G,units,MAG,PHA] = alex_tf(P,M,U)
% linearisation and numerical (Laplace) transfer function for a DCM
% witten by Alex Shaw
%
% takes a dynamical model and observation function;
%
%   dx = f(x,u,P,M)
%    y = g(x,P)
%
% and linearises flow of states [A], inputs [B] and observation [C];
%
%   dx = Ax  + Bu;
%    y = Cx [+ Du];
%
% having computed the linearisation matrices, computes the frequency response 
% at frequencies of interest (vector stored in M.Hz) using the Laplace transform.
%
% Usage: [Y,w,G,units] = alex_tf(P,M,U); 
%
% add sub-structure M.sim with M.sim.pst and M.sim.dt to force the routine
% to also recronstruct the time-domain simulation and return in the 4th
% output, since we can access the magnitude and phase data of each component
% from the Laplace transform.
%
% * Update Feb 2024: AS refactored for multiple node models to compute the
% Laplace transform of each region, then compute the cross-spectral density
%
% * Update Oct 2024: AS added second order terms; set DCM.M.tforder = 2 to
% include second order terms in transfer function (akin to including second 
% order Fourier coefficients) - although this adds significant nonlinearity
%
% AS2023

if isnumeric(P)
    P = spm_unvec(P,M.P);
end

if isstruct(P) && isfield(P,'p')
    P = P.p;
end


f = @(x,u,varargin) M.f(x,u,P,M);
w = M.Hz;

x0 = M.x(:);
u0 = 1;


% Numerically compute the Jacobian matrix
delta_x = 1e-6;  

% get delay operator
[f0,A0,D] = f(x0,u0,[]);
D         = inv(eye(length(D)) - D);

% dynamics linearisation; numerical Jacobian - dfdx
%--------------------------------------------------------------------------
[f,A,D]  = feval(M.f,M.x,0,P,M);
f  = denan(f);
A  = denan(A);

%fun = @(x) M.f(x,0,P,M);
%A=jaco(fun,M.x,ones(size(M.x))*delta_x,0,1);
%A=denan(A);

D        = inv(eye(length(D)) - D);

if isfield(M,'tforder') && M.tforder == 2
    fun = @(x) M.f(x,0,P,M);
    A2 = jaco(fun,M.x,ones(size(M.x))*delta_x,0,2);
    A2 = denan(A2);
    A2 = A2./norm(A2);
end


A        = D*A;

% input linearisation, e.g. dfdu
%--------------------------------------------------------------------------
B = spm_diff(M.f,M.x,1,P,M,2);
B = denan(B);
n = length(f);

% observation linearisation (static)
%--------------------------------------------------------------------------
C = exp(P.J);


% separate sources from here to compute sep Laplace for each unit
Ns = size(M.x,1);

% Loop each node (aka region, source, mode, column ..)
for i = 1:Ns
    win = i:Ns:(length(A));

    AA = A(win,win);
    BB = B(win);
    
    BB;% = BB + v;
    
    % we use a static observer model anyway...
    C = exp(P.J(:));
    

    % The Laplace transfer - this replaces the MATLAB built-in version using
    % 'ss' and 'bode' with a from-scratch numerical routine;
    for j = 1:length(w)
        s = exp(P.d(2)) + 1i*2*pi*w(j);
        %Jm  = AA - s*eye(length(AA));
        Jm  = s*eye(length(AA)) - AA;

        Ym  = (Jm\BB) ;%+ (Jm\x0(win));
        MG(:,j) = Ym;
        Y   = C'*Ym;
        y(j) =  Y; 
    end

    % Laplace augmented with second order term
    if isfield(M,'tforder') && M.tforder == 2
        for j = 1:length(w)
            Jm  = AA - 1i*2*pi*w(j)*eye(length(AA)) - (0.5*A2);
            Ym  = Jm\BB;
            MG(:,j) = Ym;
            Y   = C'*Ym;
            y(j) =  Y;
        end

    end

    Y = y;

    % the matlab way (not in use now)
    %G = ss(AA, BB, diag(C), 0);  % Assuming unity output matrix
    G = [];
    
    % use Bode to get Laplace transform
    %[magnitude, phase] = bode(G,w*6.2831853); % convert radians to Hz
    
    %Y = squeeze(magnitude);
    %Y = sum(Y,1);
    Y = abs(Y);

    if isfield(M,'ham') && M.ham;
        H = hamming(length(w));
        Y = Y(:).*H(:); 
    end

    MAG{i} = (MG);%magnitude;
    PHA{i} = angle(MG)*180/pi;%phase;
    
    % Laplace is pretty smooth, parameterise granularity
    H = gradient(gradient(Y));
    Y = Y - (exp(P.d(1))*3)*H;

    % inverse generalised filtering
    H = 1 ./ (1 + (w / 10).^2);

    lambda = 0.01 * exp(P.d(3));
    Gf = conj(H) ./ (abs(H).^2 + lambda);

    Y = Gf.*Y;

    % electrode scaling
    PSD(i,:) = exp(P.L(i))*(Y);

end

% compute cross spectrum from autospectra using complex cobjugate
CSD = zeros(length(w),Ns,Ns);
for i = 1:Ns
    CSD(:,i,i) = PSD(i,:);
    for j = 1:Ns
        if i ~= j
            CSD(:,i,j) = exp(P.Lc(i)) * PSD(i,:) .* conj(PSD(j,:));
            CSD(:,j,i) = CSD(:,i,j);
        end
    end
end


% global scaling / electrode gain
Y = {CSD};

units = [];

% if continuous-time simluation was requested, compute series
if isfield(M,'sim') && nargout > 3
    pst = M.sim.pst;
    dt  = M.sim.dt;
    
    % remove sim struct and recall top func
    M = rmfield(M,'sim');
    P.J(P.J==-1000)=0;
    [~,~,~,~,MAG,PHA] = atcm.fun.alex_tf(P,M,U);

    for k = 1:Ns
        
        mag = squeeze(MAG{k});
        the = squeeze(PHA{k});
    
        for i = 1:size(mag,1)
            for j = 1:size(mag,2)    
                series{k}(i,j,:) = mag(i,j) * sin(2*pi*w(j)*pst/1000 - the(i,j) );
            end
        end
    
        S{k} = squeeze(sum(series{k},2));
        S{k} = S{k} + spm_vec(M.x(k,:,:));
        LFP(k,:) = (C)'*S{k};
    end

    units.series = S;
    units.LFP    = LFP;
    units.dt     = dt;
    units.pst    = pst;
    units.A      = A;
    units.B      = B;
    units.C      = C;
    units.D      = [];
    units.xinit  = M.x(:);
    
    units.mag    = squeeze(mag);
    units.phase  = squeeze(the);
    units.freq   = w(:);

end

return;



