function [Y,w] = alex_tf(P,M,U)
% linearisation and numerical transfer function for a DCM
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
% (ss) function to get a SYS object that can be used with Bode plot to get
% the frequency response at frequencies of interest (vector stored in
% M.Hz) using essentially the Laplace transform.
%
% AS2023

f = @(x,u,varargin) M.f(x,u,P,M);
w = M.Hz;

x0 = M.x(:);
u0 = 1;

% Numerically compute the Jacobian matrix
delta_x = 1e-6;  

% get delay operator
[f0,A0,D] = f(x0,u0,[]);
D         = inv(eye(56) - D);

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

% Inputs - 
v = zeros(56,1) + 1e-3;
v(16) = v(16) + exp(P.a(1));
v(12) = v(12) + exp(P.a(2));
v(10) = v(10) + exp(P.a(3));
v(18) = v(18) + exp(P.a(4));
v(19) = v(19) + exp(P.a(5));
v(26) = v(26) + exp(P.a(6));

B = B + v;

% we use a static observer model anyway...
C = exp(P.J(:));

% Create a transfer function
%s = tf('s');
G = ss(A, B, diag(C), 0);  % Assuming unity output matrix

% use Bode to get Laplace transform
[magnitude, phase] = bode(G,w*6.2831853); % convert Hz to radians

Y = squeeze(magnitude);
Y = sum(Y,1);

% Laplace is pretty smooth, parameterise granularity
H = gradient(gradient(Y));
Y = Y - (exp(P.d(1))*3)*H;

% global scaling / electrode gain
Y = {exp(P.L(1))*abs(Y)};

