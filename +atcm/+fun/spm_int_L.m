function [y] = spm_int_L(P,M,U)
% integrates a MIMO nonlinear system using a fixed Jacobian: J(x(0))
% FORMAT [y] = spm_int_L(P,M,U)
% P   - model parameters
% M   - model structure
% U   - input structure or matrix
%
% y   - (v x l)  response y = g(x,u,P)
%__________________________________________________________________________
% Integrates the MIMO system described by
%
%    dx/dt = f(x,u,P,M)
%    y     = g(x,u,P,M)
%
% using the update scheme:
%
%    x(t + dt) = x(t) + U*dx(t)/dt
%
%            U = (expm(dt*J) - I)*inv(J)
%            J = df/dx
%
% at input times.  This integration scheme evaluates the update matrix (U)
% at the expansion point
%
%--------------------------------------------------------------------------
%
% SPM solvers or integrators
%
% spm_int_ode:  uses ode45 (or ode113) which are one and multi-step solvers
% respectively.  They can be used for any ODEs, where the Jacobian is
% unknown or difficult to compute; however, they may be slow.
%
% spm_int_J: uses an explicit Jacobian-based update scheme that preserves
% nonlinearities in the ODE: dx = (expm(dt*J) - I)*inv(J)*f.  If the
% equations of motion return J = df/dx, it will be used; otherwise it is
% evaluated numerically, using spm_diff at each time point.  This scheme is
% infallible but potentially slow, if the Jacobian is not available (calls
% spm_dx).
%
% spm_int_E: As for spm_int_J but uses the eigensystem of J(x(0)) to eschew
% matrix exponentials and inversion during the integration. It is probably
% the best compromise, if the Jacobian is not available explicitly.
%
% spm_int_B: As for spm_int_J but uses a first-order approximation to J
% based on J(x(t)) = J(x(0)) + dJdx*x(t).
%
% spm_int_L: As for spm_int_B but uses J(x(0)).
%
% spm_int_U: like spm_int_J but only evaluates J when the input changes.
% This can be useful if input changes are sparse (e.g., boxcar functions).
% It is used primarily for integrating EEG models
%
% spm_int:   Fast integrator that uses a bilinear approximation to the
% Jacobian evaluated using spm_bireduce. This routine will also allow for
% sparse sampling of the solution and delays in observing outputs. It is
% used primarily for integrating fMRI models
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_int_L.m 6270 2014-11-29 12:04:48Z karl $
 
 
% convert U to U.u if necessary
%--------------------------------------------------------------------------
if ~isstruct(U), u.u = U; U = u;   end
try, dt = U.dt;  catch, dt = 1;    end

% Initial states and inputs
%--------------------------------------------------------------------------
try
    x   = M.x;
catch
    x   = sparse(0,1);
    M.x = x;
end

try
    u   = U.u(1,:);
catch
    u   = sparse(1,M.m);
end

% add [0] states if not specified
%--------------------------------------------------------------------------
try
    f   = spm_funcheck(M.f);
catch
    f   = @(x,u,P,M) sparse(0,1);
    M.n = 0;
    M.x = sparse(0,0);
    M.f = f;
end

 
% output nonlinearity, if specified
%--------------------------------------------------------------------------
try
    g   = spm_funcheck(M.g);
catch
    g   = @(x,u,P,M) x;
    M.g = g;
end


% dx(t)/dt and Jacobian df/dx and check for delay operator
%--------------------------------------------------------------------------
D       = 1;
if nargout(f) >= 3
    [fx, dfdx,D] = f(x,u,P,M);
    
elseif nargout(f) == 2
    [fx, dfdx]   = f(x,u,P,M);
    
else
    dfdx         = spm_cat(spm_diff(f,x,u,P,M,1)); 
end
OPT.tol = 1e-6*norm((dfdx),'inf');
p       = abs(eigs(dfdx,1,'SR',OPT));
N       = ceil(max(1,dt*p*2));
n       = spm_length(x);
Q       = (spm_expm(dt*D*dfdx/N) - speye(n,n))*spm_inv(dfdx);
 
% Alex delay operators
if isfield(P,'ID')
    del = exp(P.ID).*[2 1/4 1/2 4 1/2 4 2 2]/2;
    del = repmat(del,[1 size(M.x,3)]);
    del=1./del;
    DoDelay=1;
else
    DoDelay=0;
end

if isfield(P,'ID') && size(M.x,1) > 1
    del = (spm_vec(repmat(del,[size(M.x,1) 1])))';
end

% extrinsic inter-region delays
if size(M.x,1) > 1 && isfield(P,'delays')
   
    exdel = ( round(exp(P.A{1})) | round(exp(P.A{2})) )*4;
    exdel = -exdel.*full(exp(P.delays))/1000;
    [ns,np,nk]=size(M.x);
    %exdel = kron(ones(nk,nk),kron( ~~real(P.H(:,:,1)),exdel));
    exdel = kron(ones(nk,nk),kron( ones(np,np),exdel));
    exdel = -exdel*1000;
    %del = dt*(del+exdel);
    del = dt*(del+exdel);
%     
%     Q     = (spm_expm(dt*(exdel+D)*dfdx/N) - speye(n,n))*spm_inv(dfdx);
end


% integrate
%==========================================================================
v     = spm_vec(x);
for i = 1:size(U.u,1)
    
    % input
    %----------------------------------------------------------------------
    u  = U.u(i,:);
    
    try
        % update dx = (expm(dt*J) - I)*inv(J)*f(x,u)
        %------------------------------------------------------------------
        for j = 1:N
            %v = v + Q*f(v,u,P,M);
            if DoDelay
                %v = v + del'.*(Q*f(v,u,P,M));    
                v = v + del*(Q*f(v,u,P,M));   
            else
                v = v + Q*f(v,u,P,M);
            end
        end
        
        % simpler update:
        %d = (D*dt)*v;                % delay (decay of v)
        %v = (v + d) + dt*f(v,u,P,M); % Euler scheme with delay
        
        % output - implement g(x)
        %------------------------------------------------------------------
        y(:,i) = g(v,u,P,M);

    catch
        
        % update dx = (expm(dt*J) - I)*inv(J)*f(x,u)
        %------------------------------------------------------------------
        for j = 1:N
            x = spm_vec(x) + Q*spm_vec(f(x,u,P,M));
            x = spm_unvec(x,M.x);
        end
        
        % output - implement g(x)
        %------------------------------------------------------------------
        y(:,i) = spm_vec(g(x,u,P,M));
 
        
    end
    
end
 
% transpose
%--------------------------------------------------------------------------
y      = real(y');
