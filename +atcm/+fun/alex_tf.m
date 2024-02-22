function [Y,w,G,units] = alex_tf(P,M,U)
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
% having computed the linearisation matrices, uses MATLAB's state-space
% (ss) function to get a SYS object that can be used with Bode's method to 
% get the frequency response at frequencies of interest (vector stored in
% M.Hz) using the Laplace transform.
%
% Usage: [Y,w,G,units] = alex_tf(P,M,U); 
%
% add sub-structure M.sim with M.sim.pst and M.sim.dt to force the routine
% to also recronstruct the time-domain simulation and return in the 4th
% output.
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


% Compute A matrix numerically
A = zeros(length(x0), length(x0));
for i = 1:length(x0)
    x_perturbed = x0;
    x_perturbed(i) = x_perturbed(i) + delta_x;
    dx_perturbed = f(x_perturbed, u0, []);
    A(:, i) = (dx_perturbed - f(x0, u0, [])) / delta_x;
end

A = D*A;

% Compute B matrix numerically
B = zeros(length(x0), 1);
for i = 1:length(u0)
    u_perturbed = u0;
    u_perturbed(i) = u_perturbed(i) + delta_x;
    dx_perturbed = f(x0, u_perturbed, []);
    B(:, i) = (dx_perturbed - f(x0, u0, [])) / delta_x;
end

% separate sources from here to compute sep Laplace for each unit
Ns = size(M.x,1);

% Loop each node (aka region, source, mode, column ..)
for i = 1:Ns
    win = i:Ns:(length(A));

    AA = A(win,win);
    BB = B(win);

    % Inputs - 
    v = zeros(56,1) + 1e-3;
    v(16) = v(16) + exp(P.a(1));
    v(12) = v(12) + exp(P.a(2));
    v(10) = v(10) + exp(P.a(3));
    v(18) = v(18) + exp(P.a(4));
    v(19) = v(19) + exp(P.a(5));
    v(26) = v(26) + exp(P.a(6));
    
    BB = BB + v;
    
    % we use a static observer model anyway...
    C = exp(P.J(:));
    
    % Create a transfer function
    %s = tf('s');
    G = ss(AA, BB, diag(C), 0);  % Assuming unity output matrix
    
    % use Bode to get Laplace transform
    [magnitude, phase] = bode(G,w*6.2831853); % convert radians to Hz
    
    Y = squeeze(magnitude);
    Y = sum(Y,1);

    MAG{i} = magnitude;
    PHA{i} = phase;
    
    % Laplace is pretty smooth, parameterise granularity
    H = gradient(gradient(Y));
    Y = Y - (exp(P.d(1))*3)*H;

    PSD(i,:) = exp(P.L(i))*(Y);

end

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
if isfield(M,'sim')
    pst = M.sim.pst;
    dt  = M.sim.dt;

    for k = 1:Ns
        mag = squeeze(MAG{k});
        the = squeeze(PHA{k});
    
        for i = 1:size(mag,1)
            for j = 1:size(mag,2)    
                series{k}(i,j,:) = mag(i,j) * sin(2*pi*w(j)*pst/1000 - the(i,j) );
            end
        end
    
        S{k} = squeeze(sum(series{k},2));
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
    units.phase  = squeeze(phase);
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








