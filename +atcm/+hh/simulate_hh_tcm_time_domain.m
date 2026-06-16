function out = simulate_hh_tcm_time_domain(M,P,varargin)
%SIMULATE_HH_TCM_TIME_DOMAIN Simple Euler/RK4 time-domain simulator for tc_hh_tcm.
%
%   out = simulate_hh_tcm_time_domain(M,P)
%   out = simulate_hh_tcm_time_domain(M,P,'method','rk4','T',2,'dt',1e-4)
%
% Purpose
%   A lightweight time-domain simulator for the HH-TCM model, intended for
%   exploratory checks of stability, oscillations, rebound dynamics and
%   transient responses before using the model in a DCM/spectral fitting loop.
%
% Requirements
%   M.x must be ns x 8 x 10, as returned by upgrade_tcm_to_hh.
%   M.f should normally be @tc_hh_tcm, but this function will set it if absent.
%
% Name-value options
%   'method'          'rk4' or 'euler'             default: 'rk4'
%   'T'               duration in seconds          default: 2
%   'dt'              integration step in seconds  default: 1e-4
%   'u'               constant input vector or function handle @(t) u
%                                                    default: 0
%   'source'          source to plot               default: 1
%   'populations'     populations to plot          default: 1:8
%   'plot'            true/false                   default: true
%   'clipGates'       true/false                   default: true
%   'clipConductance' true/false                   default: true
%   'verbose'         true/false                   default: true
%
% Outputs
%   out.t             Nt x 1 time vector
%   out.x             Nt x ns x 8 x 10 state trajectory
%   out.u             Nt x nu input trajectory
%   out.method        integration method
%   out.dt            integration step
%   out.popNames      population labels
%   out.stateNames    state labels
%
% Example
%   [M,P] = upgrade_tcm_to_hh(M,P);
%   M.f   = @tc_hh_tcm;
%   out   = simulate_hh_tcm_time_domain(M,P,'method','rk4','T',1,'dt',5e-5);
%
% Step input example
%   ufun = @(t) double(t > 0.1 & t < 0.3) * 0.2;
%   out  = simulate_hh_tcm_time_domain(M,P,'u',ufun,'T',1,'dt',1e-4);
%
% AS2026

% -------------------------------------------------------------------------
% Parse options
% -------------------------------------------------------------------------
ip = inputParser;
ip.addParameter('method','rk4',@(s)ischar(s)||isstring(s));
ip.addParameter('T',2,@(x)isnumeric(x)&&isscalar(x)&&x>0);
ip.addParameter('dt',1e-4,@(x)isnumeric(x)&&isscalar(x)&&x>0);
ip.addParameter('u',0);
ip.addParameter('source',1,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
ip.addParameter('populations',1:8,@(x)isnumeric(x)&&all(x>=1)&all(x<=8));
ip.addParameter('plot',true,@islogical);
ip.addParameter('clipGates',true,@islogical);
ip.addParameter('clipConductance',true,@islogical);
ip.addParameter('verbose',true,@islogical);
ip.parse(varargin{:});
opt = ip.Results;
opt.method = lower(char(opt.method));

if ~ismember(opt.method,{'euler','rk4'})
    error('method must be ''euler'' or ''rk4''.');
end

if ~isfield(M,'f') || isempty(M.f)
    M.f = @atcm.hh.tc_hh_tcm;
end

ns = size(M.x,1);
np = size(M.x,2);
nk = size(M.x,3);

if np ~= 8 || nk < 10
    error('Expected M.x to be ns x 8 x 10. Run [M,P] = upgrade_tcm_to_hh(M,P) first.');
end
if opt.source > ns
    error('Requested source %d but M.x only has %d source(s).',opt.source,ns);
end

Nt = floor(opt.T/opt.dt) + 1;
t  = (0:Nt-1)' * opt.dt;

xvec = spm_vec(M.x);
Nx   = numel(xvec);

% Determine input dimension from first u call.
u0 = local_get_u(opt.u,t(1));
u0 = u0(:);
Nu = numel(u0);

X = zeros(Nt,Nx);
U = zeros(Nt,Nu);
X(1,:) = xvec(:)';
U(1,:) = u0(:)';

if opt.verbose
    fprintf('HH-TCM simulation: %s, T = %.4g s, dt = %.4g s, steps = %d, states = %d\n', ...
        opt.method,opt.T,opt.dt,Nt,Nx);
end

% -------------------------------------------------------------------------
% Main integration loop
% -------------------------------------------------------------------------
for n = 1:Nt-1
    tn = t(n);
    un = local_get_u(opt.u,tn);
    un = un(:);
    U(n,:) = un(:)';

    switch opt.method
        case 'euler'
            dx = local_rhs(xvec,un,P,M);
            xvec = xvec + opt.dt * dx;

        case 'rk4'
            u1 = local_get_u(opt.u,tn);                  u1 = u1(:);
            u2 = local_get_u(opt.u,tn + opt.dt/2);       u2 = u2(:);
            u3 = u2;
            u4 = local_get_u(opt.u,tn + opt.dt);         u4 = u4(:);

            k1 = local_rhs(xvec,                 u1,P,M);
            k2 = local_rhs(xvec + opt.dt*k1/2,   u2,P,M);
            k3 = local_rhs(xvec + opt.dt*k2/2,   u3,P,M);
            k4 = local_rhs(xvec + opt.dt*k3,     u4,P,M);
            xvec = xvec + opt.dt*(k1 + 2*k2 + 2*k3 + k4)/6;
    end

    xvec = local_sanitise_state(xvec,M,opt.clipGates,opt.clipConductance);

    if any(~isfinite(xvec))
        warning('Non-finite state encountered at step %d, t = %.6f s. Truncating output.',n,tn);
        X = X(1:n,:);
        U = U(1:n,:);
        t = t(1:n);
        Nt = n;
        break
    end

    X(n+1,:) = xvec(:)';
end

U(end,:) = local_get_u(opt.u,t(end));

% Reshape to Nt x ns x np x nk
X4 = reshape(X,[Nt ns np nk]);

out = struct();
out.t = t;
out.x = X4;
out.u = U;
out.method = opt.method;
out.dt = opt.dt;
out.M = M;
out.P = P;
out.popNames = {'SS','SP','SI','DP','DI','TP','RT','RC'};
out.stateNames = {'V','gA','gG','gN','gB','pM','rH','mNa','hNa','nK'};

if opt.plot
    local_plot_out(out,opt.source,opt.populations);
end
end

% =========================================================================
% Local functions
% =========================================================================
function dx = local_rhs(xvec,u,P,M)
dx = M.f(xvec,u,P,M);
dx = dx(:);
end

function u = local_get_u(uSpec,t)
if isa(uSpec,'function_handle')
    u = uSpec(t);
else
    u = uSpec;
end
if isempty(u)
    u = 0;
end
end

function xvec = local_sanitise_state(xvec,M,clipGates,clipConductance)
ns = size(M.x,1);
np = size(M.x,2);
nk = size(M.x,3);
x  = reshape(xvec,ns,np,nk);

% Conductances should not become negative during exploratory simulation.
if clipConductance
    x(:,:,2:5) = max(x(:,:,2:5),0);
end

% Gates are probabilities/activation variables.
if clipGates && nk >= 10
    x(:,:,6:10) = min(max(x(:,:,6:10),0),1);
end

xvec = x(:);
end

function local_plot_out(out,src,pops)
t = out.t;
X = out.x;
popNames = out.popNames;

V    = squeeze(X(:,src,:,1));
gA   = squeeze(X(:,src,:,2));
gG   = squeeze(X(:,src,:,3));
gN   = squeeze(X(:,src,:,4));
gB   = squeeze(X(:,src,:,5));
pM   = squeeze(X(:,src,:,6));
rH   = squeeze(X(:,src,:,7));
mNa  = squeeze(X(:,src,:,8));
hNa  = squeeze(X(:,src,:,9));
nK   = squeeze(X(:,src,:,10));

% Simple firing proxy matching tc_hh_tcm default, using VR0 + 0 if P.S is not
% easily recoverable from the output structure.
VR0 = -52;
R   = 0.25;
mfireApprox = 1 ./ (1 + exp(-R .* (V - VR0)));

figure('Name','HH-TCM time-domain voltage','Color','w');
plot(t,V(:,pops),'LineWidth',1.1);
xlabel('Time (s)'); ylabel('Voltage (mV)');
title(sprintf('HH-TCM voltage, source %d',src));
legend(popNames(pops),'Location','bestoutside'); box off;

figure('Name','HH-TCM approximate firing proxy','Color','w');
plot(t,mfireApprox(:,pops),'LineWidth',1.1);
xlabel('Time (s)'); ylabel('Approx. firing proxy');
title(sprintf('HH-TCM firing proxy, source %d',src));
legend(popNames(pops),'Location','bestoutside'); box off;

figure('Name','HH-TCM synaptic conductances','Color','w');
plot(t,gA(:,pops),'LineWidth',1.0); hold on;
plot(t,gG(:,pops),'LineWidth',1.0);
plot(t,gN(:,pops),'LineWidth',1.0);
plot(t,gB(:,pops),'LineWidth',1.0);
xlabel('Time (s)'); ylabel('Conductance state');
title(sprintf('Synaptic conductance states, source %d',src));
box off;

figure('Name','HH-TCM active gates','Color','w');
plot(t,mNa(:,pops),'LineWidth',1.0); hold on;
plot(t,hNa(:,pops),'LineWidth',1.0);
plot(t,nK(:,pops),'LineWidth',1.0);
plot(t,pM(:,pops),'LineWidth',1.0);
plot(t,rH(:,pops),'LineWidth',1.0);
xlabel('Time (s)'); ylabel('Gate value');
title(sprintf('HH/M/H gates, source %d',src));
box off;

figure('Name','HH-TCM thalamic focus','Color','w');
thpops = intersect(pops,[6 7 8]);
if isempty(thpops), thpops = [6 7 8]; end
plot(t,V(:,thpops),'LineWidth',1.2);
xlabel('Time (s)'); ylabel('Voltage (mV)');
title(sprintf('Thalamo-cortical populations, source %d',src));
legend(popNames(thpops),'Location','bestoutside'); box off;
end
