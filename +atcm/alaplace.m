function [y,w,yt,t,ys] = alaplace(P,M,U,varargin)
% Numerical Laplace transform for TCM. This version produces the same result 
% as 'atcm.fun.alex_tf', except this version implements tghe Laplace transform 
% from scratch instead of relying on the matlab functions ss and bode.
%
% AS2024

if nargin == 4 && ~isempty(varargin{1})
    %fprintf('using user-supplied M.x\m');
    M.x = spm_unvec(varargin{1},M.x);
    x = M.x;
end

w = M.Hz;

% solve for steady-state - if required
%--------------------------------------------------------------------------   
%M.x   = atcm.fun.solvefixedpoint(P,M,0,-70);


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


% Laplace transform
%--------------------------------------------------------------------------
for j = 1:length(w)

    Jm = A - 1i*2*pi*w(j)*eye(n);
    Y  = Jm\B;
    %Y  = M.x(:)'*Y;
    Y = C'*Y;

    y(j) =  Y; 

end

% Laplace is pretty smooth, parameterise granularity
H = gradient(gradient(y));
y = y - (exp(P.d(1))*3)*H;

y = exp(P.L) * y;
y = {abs(y(:))};

% numerical integration
if nargout < 3; return; end


t  = M.sim.pst;
dt = M.sim.dt;

v       = M.x(:);
yt(:,1) = v;

for j = 1:length(t)

    dxdt = feval(M.f,v,0,P,M) ;
    b    = pinv(full(A)'.*v).*dxdt;
    b    = b + B'*dt;
    Q    = A.*b; 
    v    = v + dt*Q*v ;

    yt(:,j) = v;
end

% S  = C'*yt;
 fs = 1/dt;

% 
% 
% %Fs = atcm.fun.asinespectrum(w,t,[],@sin);
% %Fs = Fs + 1i*2*pi * Fs;
% %X   = [ones(1,size(Fs,2)); Fs];
% 
for i = 1:size(yt,1)
    Pf(i,:) = atcm.fun.Afft(yt(i,:),fs,w);
end 

ys  = abs(C'*Pf);

ys = atcm.fun.agauss_smooth(ys,1);

% Laplace is pretty smooth, parameterise granularity
%H = gradient(gradient(ys));
%ys = ys - (exp(P.d(1))*3)*H;
ys = exp(P.L) * ys;


