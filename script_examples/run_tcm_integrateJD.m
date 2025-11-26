function [X, time, info, u_t] = run_tcm_integrateJD(M, P, varargin)
%RUN_TCM  Convenience wrapper to integrate a TCM/DCM with static-JD integrator.
%
%   [X, t, info, u_t] = run_tcm(M, P, Name,Value, ...)
%
% Required
%   M    : model struct (expects fields .f (handle) and .x (initial states))
%   P    : parameter struct (passed through to M.f and integrator)
%
% Name-Value options (all optional)
%   'dt'          : time step (s). Default 1/600.
%   'T'           : total duration (s). Default 2.
%   't'           : explicit time vector; overrides dt/T if provided.
%   'x0'          : initial state vector. Default M.x(:) if present, else zeros.
%
%   Input options (choose one path below):
%   'u'           : custom input. Either function handle @(tt)->scalar/vec,
%                   or numeric vector (same length as t). Overrides other input opts.
%   'amp'         : amplitude for the default sine burst. Default 1/128.
%   'freq'        : frequency (Hz) for sine burst. Default 20.
%   't0'          : burst onset time (s). Default 0.
%   'dur'         : burst duration (s). Default 0.5.
%
%   Integrator options (passed to atcm.tcm_integrate_staticJD):
%   'method'      : integrator method. Default 'rosenbrock'.
%   'substeps'    : integer sub-steps. Default 8.
%   'autoAlpha'   : logical. Default true.
%   'finiteCheck' : logical. Default true.
%
% Returns
%   X    : state trajectories (nStates × nTime)
%   t    : time vector (1 × nTime)
%   info : struct returned by the integrator
%   u_t  : input evaluated on t (1 × nTime), if scalar input; [] otherwise
%
% Example
%   [X,t,info,u] = run_tcm(DCM.M, DCM.M.pE, 'dt',1/600, 'T',2, 'freq',20, 'dur',0.5);
%

% x = reshape(X,[1 8 7 1200]);lab = {'mV' 'AMPA' 'GABA-A' 'NMDA' 'GABA-B' 'Kv7 (M)' 'HCN'};
% for i = 1:7
%     subplot(7,1,i);
%     plot(t,squeeze(x(1,:,i,:)));title(lab{i});
% end

% -------------------- parse inputs --------------------
ip = inputParser;
ip.FunctionName = 'run_tcm_integrateJD';

addParameter(ip, 'dt',          1/600, @(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(ip, 'T',           2.0,   @(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(ip, 'time',           [],    @(x)isnumeric(x)&&isvector(x));
addParameter(ip, 'x0',          [],    @(x)isnumeric(x)&&isvector(x));
addParameter(ip, 'u',           [],    @(x)isa(x,'function_handle')||isnumeric(x)||isempty(x));

addParameter(ip, 'amp',         1/128, @(x)isnumeric(x)&&isscalar(x));
addParameter(ip, 'freq',        20,    @(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(ip, 't0',          0.0,   @(x)isnumeric(x)&&isscalar(x));
addParameter(ip, 'dur',         0.5,   @(x)isnumeric(x)&&isscalar(x)&&x>=0);

addParameter(ip, 'method',      'rosenbrock', @(s)ischar(s)||isstring(s));
addParameter(ip, 'substeps',    8,     @(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(ip, 'autoAlpha',   true,  @(x)islogical(x)||ismember(x,[0 1]));
addParameter(ip, 'finiteCheck', true,  @(x)islogical(x)||ismember(x,[0 1]));

parse(ip, varargin{:});
prm = ip.Results;

% -------------------- time vector --------------------
if isempty(prm.time)
    time = 0:prm.dt:(prm.T - prm.dt);
else
    time = prm.time(:)';  % ensure row
end

% -------------------- initial state --------------------
if isempty(prm.x0)
    if isfield(M,'x') && ~isempty(M.x)
        x0 = M.x(:);
    else
        % fallback: try M.n or infer from f at zeros
        if isfield(M,'n') && ~isempty(M.n)
            x0 = zeros(M.n,1);
        else
            error('run_tcm:NoInitialState', ...
                ['Initial state not provided and M.x/M.n not available.\n' ...
                 'Pass ''x0'', or ensure M.x exists.']);
        end
    end
else
    x0 = prm.x0(:);
end

% -------------------- build input u(tt) --------------------
u_t = [];  % default return for input sampled on grid (if scalar)
if ~isempty(prm.u)
    % custom input provided
    if isa(prm.u, 'function_handle')
        ufun = prm.u;
    else
        % numeric vector: make time-locked piecewise-constant/linear interpolant
        uvec = prm.u(:)';  % row
        if numel(uvec) ~= numel(time)
            error('run_tcm:BadU', 'Numeric u must have same length as t (%d).', numel(time));
        end
        % piecewise-constant (fast & typical for inputs); change to 'linear' if preferred
        ufun = @(tt) interp1(time, uvec, tt, 'previous', 'extrap');
        u_t  = uvec;
    end
else
    % default: windowed sine burst
    amp  = prm.amp;
    f    = prm.freq;
    t0   = prm.t0;
    dur  = prm.dur;
    ufun = @(tt) double(tt >= t0 & tt < (t0 + dur)) .* (amp .* sin(2*pi*f*(tt - t0)));
    % also sample it on the grid (assuming scalar input)
    u_t  = ufun(time);
end

% -------------------- integrator options --------------------
opts = struct( ...
    'method',      char(prm.method), ...
    'substeps',    prm.substeps, ...
    'autoAlpha',   logical(prm.autoAlpha), ...
    'finiteCheck', logical(prm.finiteCheck) ...
    );

% -------------------- call integrator --------------------
if ~isfield(M,'f') || isempty(M.f) || ~isa(M.f,'function_handle')
    error('run_tcm:NoF', 'M.f must be a valid function handle: dx = f(x,u,P,M)');
end

[X, info] = atcm.tcm_integrate_staticJD(M.f, M, P, x0, time, ufun, opts);
end
