function DCM = dcm_erp_invert(DCM)

M  = DCM.M;
xU = DCM.xU;
xY = DCM.xY;

if strcmp( char(DCM.M.IS), 'atcm.integrate_erp')
    M.IS = 'atcm.integrate_erp'; % remove '@'
end
fprintf('Integrator: %s\n',char(DCM.M.IS));
fprintf('Model function: %s\n',char(DCM.M.f));

M.nograph = 0;

% EM: inversion
%==========================================================================
[Qp,Qg,Cp,Cg,Ce,F,LE] = spm_nlsi_N(M,xU,xY);

% Data ID
%==========================================================================
if isfield(M,'FS')
    try
        ID = spm_data_id(feval(M.FS,xY.y,M));
    catch
        ID = spm_data_id(feval(M.FS,xY.y));
    end
else
    ID = spm_data_id(xY.y);
end

pE = DCM.M.pE;
pC = DCM.M.pC;

% Bayesian inference
%--------------------------------------------------------------------------
sw  = warning('off','SPM:negativeVariance');
dp  = spm_vec(Qp) - spm_vec(pE);
Pp  = spm_unvec(1 - spm_Ncdf(0,abs(dp),diag(Cp)),Qp);
warning(sw);


% neuronal and sensor responses (x and y)
%--------------------------------------------------------------------------
L   = feval(M.G, Qg,M);                 % get gain matrix
x   = feval(M.IS,Qp,M,xU);              % prediction (source space)


% trial-specific responses (in mode, channel and source space)
%--------------------------------------------------------------------------
Ns  = length(DCM.A{1}); 
Nr  = size(DCM.C,1); 
Nt  = length(x);
j   = find(kron( exp(Qg.J),               ones(1,Nr)));      % Indices of contributing states
%j    = find(kron( logical(sum(Qg.J)), ones(1,Nr)));

NNs = size(xY.y{1},1);
x0  = ones(NNs,1)*spm_vec(M.x)';         % expansion point for states
for i = 1:Nt
    K{i} = x{i} - x0;                   % centre on expansion point
    y{i} = M.R*K{i}*L'*M.U;             % prediction (sensor space)
    r{i} = M.R*xY.y{i}*M.U - y{i};      % residuals  (sensor space)
    K{i} = K{i}(:,j);                   % Depolarization in sources
end


% store estimates in DCM
%--------------------------------------------------------------------------
DCM.M  = M;                    % model specification
DCM.xY = xY;                   % data structure
DCM.xU = xU;                   % input structure
DCM.Ep = Qp;                   % conditional expectation f(x,u,p)
DCM.Cp = Cp;                   % conditional covariances G(g)
DCM.Eg = Qg;                   % conditional expectation
DCM.Cg = Cg;                   % conditional covariances
DCM.Ce = Ce;                   % conditional error
DCM.Pp = Pp;                   % conditional probability
DCM.H  = y;                    % conditional responses (y), projected space
DCM.K  = K;                    % conditional responses (x) (contributing)
DCM.x  = x;                    % conditional responses (x) (all states)
DCM.R  = r;                    % conditional residuals (y)
DCM.F  = F;                    % Laplace log evidence
DCM.L  = LE;                   % Laplace log evidence components
DCM.ID = ID;                   % data ID

DCM.options.Nmodes   = size(M.U,2);
%DCM.options.onset    = onset;
%DCM.options.dur      = dur;
%DCM.options.model    = model;
%DCM.options.lock     = lock;
%DCM.options.symm     = symm;
DCM.options.analysis = 'ERP';


% and save
%--------------------------------------------------------------------------
[fp fn fe] = fileparts(DCM.name);
DCM.name = [fn fe];
save(DCM.name, 'DCM', spm_get_defaults('mat.format'));


end