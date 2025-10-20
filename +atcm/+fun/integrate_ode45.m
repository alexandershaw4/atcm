function [y,t] = integrate_ode45(P,M,U,pst)
% integrate a neuronal model of DCM-like form using a matlab built-in integrator
%
%
%
% AS

% expansion (fixed) point: trial & parameter effects are deviations from here
%--------------------------------------------------------------------------
f    = spm_funcheck(M.f); 
x    = atcm.fun.solvefixedpoint(P,M);
M.x  = x;
v    = spm_vec(x);
NoFX = 0;

if isempty(U)
    U.X  = 1;
    NoFX = 1; % flag no modulatory effects in this model
end

dt    = M.dt;
Fs    = 1/dt;                     % sampling frequency
%tn    = 2;                        % sample window length, in seconds
%pst   = 1000*((0:dt:tn-dt)');     % peristim times we'll sample
tspan = pst;
ns    = length(pst);

U.dt = dt;

% output nonlinearity
%--------------------------------------------------------------------------
% try
%     g   = fcnchk(M.g,'x','u','P','M');
% catch
%     g   = [];
% end

for  c = 1:size(U.X,1)
    
    % generate condition-specific parameter structure
    %----------------------------------------------------------------------
    if ~NoFX; Q  = spm_gen_Q(P,U.X(c,:));
    else      Q  = P;
    end
        
    % integration, spectral response, firing etc. (subfunction)
%     ode   = inline('spm_vec(f(spm_unvec(x,M.x),0,Q,M))',...
%                't','x','OPTIONS','Q','M','U','f');
    
    % Newer version require an anonymous function
    ode = @(t,x,Q,M,f) spm_vec( f(spm_unvec(x,M.x),0,Q,M) );

    %OPTIONS = odeset;
    
    % Runge-Kutta 45 scheme
    [t,x]   = ode113(ode,tspan/1000,spm_vec(x),' ',Q,M,f);
    
    % Compute Lyapunov Exponents
    %[Texp,Lexp]=lyapunov(56,ode,@ode113,1,dt,2,spm_vec(x),1,{' ',Q,M,f});
    
%     % output
%     %----------------------------------------------------------------------
%     for i = 1:ns
% 
%         % output - implement g(x)
%         %------------------------------------------------------------------
%         if ~isempty(g)
%             y(:,i) = spm_vec(g(spm_unvec(x(i,:),M.x),U.u(i,:),P,M));
%         else
%             y(:,i) = spm_vec(spm_unvec(x(i,:),M.x));
%         end
% 
%     end
% 
%     % transpose
%     %----------------------------------------------------------------------
%     y      = real(y');
    
y = x;
    
    
end



           


