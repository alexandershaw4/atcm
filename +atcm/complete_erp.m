function DCM = complete_erp(DCM)

DCM.CUSTOM.nograph   = 0;
DCM.options.multiC   = 0;

% Filename and options
%--------------------------------------------------------------------------
try, Nm       = DCM.options.Nmodes;   catch, Nm        = 8;           end
try, onset    = DCM.options.onset;    catch, onset     = 60;          end
try, dur      = DCM.options.dur;      catch, dur       = 16;          end
try, model    = DCM.options.model;    catch, model     = 'NMM';       end
try, lock     = DCM.options.lock;     catch, lock      = 0;           end
try, multC    = DCM.options.multiC;   catch, multC     = 0;           end
try, symm     = DCM.options.symmetry; catch, symm      = 0;           end
try, CVA      = DCM.options.CVA;      catch, CVA       = 0;           end
try, Nmax     = DCM.options.Nmax;     catch, Nmax      = 64;          end
try, DATA     = DCM.options.DATA;     catch, DATA      = 1;           end
try, Nmax     = DCM.M.Nmax;           catch, Nmax      = Nmax;        end


% symmetry contraints for ECD models only
%--------------------------------------------------------------------------
if ~strcmp(DCM.options.spatial,'ECD'), symm = 0; end


% Data and spatial model
%==========================================================================
%DCM = spm_dcm_erp_data(DCM);
DCM = spm_dcm_erp_dipfit(DCM,1);
xY  = DCM.xY;
xU  = DCM.xU;
M   = DCM.M;

% dimensions
%--------------------------------------------------------------------------
Nt      = length(xY.y);                  % number of trials
Nr      = size(DCM.C,1);                 % number of sources
Nu      = size(DCM.C,2);                 % number of exogenous inputs
Ns      = size(xY.y{1},1);               % number of time bins
Nc      = size(xY.y{1},2);               % number of channels
Nx      = size(xU.X,2);                  % number of trial-specific effects

% check the number of modes is greater or equal to the number of sources
%--------------------------------------------------------------------------
Nm      = max(Nm,Nr);

% confounds - residual forming matrix
%--------------------------------------------------------------------------
if isfield(xY,'R')
    M.R = xY.R;
else
    M.R = speye(Ns) - xY.X0*((xY.X0'*xY.X0)\xY.X0');
end


% Serial correlations (precision components) AR model
%--------------------------------------------------------------------------
xY.Q   = {spm_Q(1/2,Ns,1)};

% between-trial effects
%--------------------------------------------------------------------------
try
    if length(DCM.B) < Nx
        for i = 1:Nx
            DCM.B{i} = sparse(Nr,Nr);
        end
    end
catch
    xU.X  = sparse(1,0);
    DCM.B = {};
end

% within-trial effects: adjust onset relative to PST
%--------------------------------------------------------------------------
M.ons  = onset - xY.pst(1);
M.dur  = dur;
xU.dt  = xY.dt;


%-Model specification and nonlinear system identification
%==========================================================================
try, M = rmfield(M,'g'); end

% prior moments on parameters
%--------------------------------------------------------------------------
pE = DCM.M.pE;
pC = DCM.M.pC;

% priors on spatial model
%--------------------------------------------------------------------------
M.dipfit.model = model;
%[gE,gC]        = spm_L_priors(M.dipfit);


% Copy the (spatial) priors into their own structure
gE.J    = pE.J;
gE.L    = pE.L;
gE.Lpos = pE.Lpos;

gC.J    = pC.J;
gC.L    = pC.L;
gC.Lpos = pC.Lpos;

% Set prior correlations (locking trial effects and dipole orientations
%--------------------------------------------------------------------------
if lock, pC = spm_dcm_lock(pC);      end
if symm, gC = spm_dcm_symm(gC,gE);   end


% hyperpriors (assuming a high signal to noise)
%--------------------------------------------------------------------------
hE      = 6;
hC      = 1/128;

% scale data features
%--------------------------------------------------------------------------
scale    = std(spm_vec(spm_fy_erp(xY.y,M)));
xY.y     = spm_unvec(spm_vec(xY.y)/scale,xY.y);
xY.scale = xY.scale/scale;

% Handles
IS = DCM.M.IS;

fprintf('\nInitialising model (hidden) states\n');
%[x,f,h] = spm_dcm_x_neural(pE,model);

Q   = atcm.fun.spm_gen_Q_as(pE,DCM.xU.X(1,:));
[x] = atcm.fun.solvefixedpoint(Q,M);
h=1;

M.FS   = 'spm_fy_erp';
%M.G    = 'spm_lx_erp_tc6';
%M.G = @atcm.fun.spm_lx_erp_tc6;
M.G = 'spm_lx_erp';
M.IS   = IS;
M.h    = h;
M.x    = x;
M.pE   = pE;
M.pC   = pC;
M.gE   = gE;
M.gC   = gC;
M.hE   = hE;
M.hC   = hC;
M.m    = Nu;
M.n    = length(spm_vec(M.x));
M.l    = Nc;
M.ns   = Ns;
M.Nmax = Nmax;

% Complete specification for use with aINTw_net integration
M.Hz     = 4:80;
M.pst    = xY.pst / 1000;
M.dt     = xY.dt;

% overwrite neural priors - for now (until we have better)
%fprintf('Re-setting priors\n');
%M.pE = spm_unvec( spm_vec(M.pE)*0, M.pE);

% Check that f and G are callable, from the equation:
%  { y  = G(dx,P,M)      = spatial projection
%  { dx = IS(f(x,u,P,M)) = numerical integration of neural model f
fprintf('Checking integration and spatial projection functions\n');

L   = feval(M.G, gE,M);            
x   = feval(M.IS,pE,M,xU);  

% return everything
DCM.M  = M;
DCM.xU = xU;
DCM.xY = xY;

% To reconstruct the model (initial) predicted timeseries for trial j;
% Data = L'*x{j}
%
% & plot:
% plot(xY.pst,Data)

