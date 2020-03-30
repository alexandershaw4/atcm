function DCM = complete(DCM)
% Complete specification of the tcm model ready for inversion
%
%
%

DCM.CUSTOM.nograph   = 0;
DCM.options.multiC   = 0;

% Filename and options
%--------------------------------------------------------------------------
try, model   = DCM.options.model;   catch, model    = 'NMM';     end
try, spatial = DCM.options.spatial; catch, spatial  = 'LFP';     end
try, Nm      = DCM.options.Nmodes;  catch, Nm       = 8;         end
try, DATA    = DCM.options.DATA;    catch, DATA     = 1;         end


% Spatial model
%==========================================================================
DCM.options.Nmodes = Nm;
DCM.M.dipfit.model = model;
DCM.M.dipfit.type  = spatial;

if ~DATA 
    DCM  = spm_dcm_erp_data(DCM);                   % data
end
if DATA
    DCM  = spm_dcm_erp_dipfit(DCM, 1);              % spatial model
end
Ns   = length(DCM.A{1});                            % number of sources


% Design model and exogenous inputs
%==========================================================================
if ~isfield(DCM,'xU'),   DCM.xU.X = sparse(1 ,0); end
if ~isfield(DCM.xU,'X'), DCM.xU.X = sparse(1 ,0); end
if ~isfield(DCM,'C'),    DCM.C    = sparse(Ns,0); end
if isempty(DCM.xU.X),    DCM.xU.X = sparse(1 ,0); end
if isempty(DCM.xU.X),    DCM.C    = sparse(Ns,0); end

% Neural mass model parameters
%==========================================================================
pE = DCM.M.pE;

% initial states and equations of motion
%--------------------------------------------------------------------------
[x,~] = atcm.fun.solvefixedpoint(pE,DCM.M);
 
% create DCM
%--------------------------------------------------------------------------
DCM.M.FS = 'spm_fs_csd';
DCM.M.g  = 'spm_gx_erp';
%DCM.M.f  = f;
DCM.M.x  = x;
DCM.M.n  = length(spm_vec(x));
DCM.M.hE = 6;
DCM.M.hC = 1/64;
DCM.M.m  = Ns;

% specify M.u - endogenous input (fluctuations) and intial states
%--------------------------------------------------------------------------
DCM.M.u  = sparse(Ns,1);

%-Feature selection using principal components (U) of lead-field
%==========================================================================
 % re-initialise states given new priors
[x,f] = atcm.fun.solvefixedpoint(pE,DCM.M);
DCM.M.x = x;

% get data-features (in reduced eigenspace)
%==========================================================================
if ~DATA
   % DCM  = spm_dcm_csd_data(DCM);
   try fdata = DCM.options.data;
   catch
       fdata = 'spm_dcm_csd_data_as';
   end
   fdata
    %DCM  = spm_dcm_csd_data_as(DCM);
    DCM = feval(fdata,DCM);
end
 
% scale data features (to a variance of about 8)
%--------------------------------------------------------------------------
ccf      = spm_csd2ccf(DCM.xY.y,DCM.xY.Hz);
scale    = max(spm_vec(ccf));
DCM.xY.y = spm_unvec(8*spm_vec(DCM.xY.y)/scale,DCM.xY.y);



% complete model specification and invert
%==========================================================================
Nm       = size(DCM.M.U,2);                    % number of spatial modes
DCM.M.l  = Nm;
DCM.M.Hz = DCM.xY.Hz;
DCM.M.dt = DCM.xY.dt;
 
% precision of noise: AR(1/2)
%--------------------------------------------------------------------------
% y     = spm_fs_csd(DCM.xY.y,DCM.M);
% DCM.M.FS = 'spm_fs_csd';

% Alex new features selection - just abs(csd)
%--------------------------------------------------------------------------
FS = inline('spm_unvec(real(spm_vec(x)),x)','x','y');
y  = feval(FS,DCM.xY.y,DCM.M);
DCM.M.FS = FS;

for i = 1:length(y)
    n      = size(y{i},1);
    m      = size(y{i},2)*size(y{i},3);
    q      = spm_Q(1/2,n,1);  
    
    % alex - focus on higher frequencies: Q*diag(Hz)*Q
    if  n == length(DCM.M.Hz)
        H = DCM.M.Hz;
    else
        H = linspace(DCM.M.Hz(1),DCM.M.Hz(end),n);
    end
    
    q = spm_Q(1/2,n,1)*diag(H)*spm_Q(1/2,n,1); % as per ssr gamma
    Q{i,i} = kron(speye(m,m),q);
end
DCM.xY.Q  = spm_cat(Q);
DCM.xY.X0 = sparse(size(Q,1),0);



try NGP = DCM.CUSTOM.nograph;
    if NGP
        DCM.M.nograph = 1;
    end
end
try DCM.M.ShowNodes = CUS.ShowNodes; end