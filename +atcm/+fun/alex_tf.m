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
% AS2023

if isnumeric(P)
    P = spm_unvec(P,M.P);
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
J        = A;
D        = inv(eye(length(D)) - D);
A        = D*A;

% input linearisation, e.g. dfdu
%--------------------------------------------------------------------------
B = spm_diff(M.f,M.x,1,P,M,2);
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
        Jm  = AA - 1i*2*pi*w(j)*eye(length(AA));
        Ym  = Jm\BB;
        MG(:,j) = Ym;
        Y   = C'*Ym;
        y(j) =  Y;
    end

    Y = y;

    % the matlab way (not in use now)
    G = ss(AA, BB, diag(C), 0);  % Assuming unity output matrix
    
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

    PSD(i,:) = exp(P.L(i))*(Y);

end

% compute cross spectrum from autospectra using complex cobjugate
CSD = zeros(length(w),Ns,Ns);
for i = 1:Ns
    CSD(:,i,i) = PSD(i,:);
    for j = 1:Ns
        if i ~= j
            CSD(:,i,j) = PSD(i,:) .* conj(PSD(j,:));
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
        LFP(k,:) = diag(G.C)'*S{k};
    end

    units.series = S;
    units.LFP    = LFP;
    units.dt     = dt;
    units.pst    = pst;
    units.A      = G.A;
    units.B      = G.B;
    units.C      = G.C;
    units.D      = G.D;
    units.xinit  = M.x(:);
    
    units.mag    = squeeze(mag);
    units.phase  = squeeze(the);
    units.freq   = w(:);

end

return;

% % notes on how to use this model to integrate in the time domain but remove
% % log-linear trend;
% %--------------------------------------------------------------------------
% 
A = G.A;
B = G.B;
C = G.C;

dt = 1/1000;

% initial point;
x = dt*A*M.x(:) + dt*B;

% Euler
for i = 2:100; 
    x(:,i) = x(:,i-1) + dt*A*x(:,i-1) + B; 
end
% 
% % log-linear trend
% for i = 1:size(x,1)
%     proj = log(x(i,:));
%     Y = fit(t.',proj.','poly1');
%     xx(i,:) = exp( proj - Y(t)' );
% end
% 
% Yactual = diag(C)'*xx;








