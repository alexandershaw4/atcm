function [X0,X,cm,F,j] = icaoptimse(P,DCM,method,nc)
global DD
% Optimisation / inversion of DCM models using parameter space partitioning
% and optimisation via minimising either:
%   e = sum( data - model ).^2 
%
% or the free-energy.
%   see atcm.optim.fe0.m
%
% Parameter space paritioning (reduction) is achieved by 
% clustering of the system Jacobian (df/dx) for parameters of non-zero
% prior variance. This requires a one-off calculation of J regardless of
% subsequent optim method. Set nc == the number of non-zero priors to skip
% this step. This partitioning is conceptually a 'surrogate model'
% approach. The subsequent optimisation (of this surrogate model) may be 
% either gradient based or gradient free.
%
% Optimisation is achieved either using the default gradient-based DCM Bayesian
% EM routine (spm_nlsi_GN), or through minimisation of the objective function using:
%
% - fminsearch (gradient based)
% - hooke-jeeves pattern search algorithm (gradient free)
% - simulated annealing (gradient free, requires global optimisation toolbox)
% - genetic algorithm (requires global optimisation toolbox)
% - particle swarm
% - patternsearch
% - laplacian sampling (see spm_nlsi_LS)
% - surrogate model optimisation (requires global optimisation toolbox)
%
% Function usage:
%  [X0,X,cm,F,j] = clusteroptimse(P,DCM,method,nc)
%
% P   = parameter set in DCM structure, e.g. DCM.M.pE or DCM.Ep
% DCM = the main DCM data structure, containing:
%       - DCM.M.IS = name / handle of the integration function 
%       - DCM.M.f  = name / handle of the equations of motion function
%       - DCM.M.pC = prior variances structure corresponding to P
%
% method = switch for optim method, options:
%     'abc'       - artificial bee colony optimsation
%     'fminsearch' - built in nonlinear solver for ~ 6 param problems
%     'hj'        - hooke-jeeves pattern search 
%     'sa'        - simulated annealing for [overspecified problems]
%     'ga'        - genetic algorithm
%     'nlsi'      - Default DCM Bayesian Expentation-Maximisation routine
%     'particle'  - particle swarm optimsations using global optim toolbox
%     'pattern'   - generalised apptern search using global optim toolbox
%     'ls'        - DCM Bayesian EM algorithm but with Laplacian sampling
%     'surrogate' - surrogate model optimisation; attempts to build a
%     surrogate function approximating the full system
%
% nc     = number of components (params) to reduce the full param set to.
% (nc can be a vector to iterate over)
%
% AS2018

DD   = DCM;
DD.P = P;

% number of iterations (re-calc of reduced space, re-fit)
%--------------------------------------------------------------------------
nr = length(nc);    

% feature selection: params of max effect, or hierarchical clustering 
%--------------------------------------------------------------------------
fs = @clusteroptimse0; %{'clusteroptimse_maxi' or 'clusteroptimse0'};

% include (prior) estimation of parameter covariance (or not)
%--------------------------------------------------------------------------
pC = [];      % do not use prior covariance estimate
%pC = DCM.Cp; % parameter covariance matrix


for i = 1:nr
    fprintf('Iteration %d / %d (with %i components)\n',i,nr,nc(i));
    
    % Compute a parameter set in reduced space
    %----------------------------------------------------------------------
    [Px, cm{i}, cx{i}, j{i}] = fs(P,DCM,nc(i),pC);
    DD.cm                    = cm{i};
    DD.cx                    = cx{i};

    
    % Optimise the model.... options:
    %----------------------------------------------------------------------
    switch method
        
        case 'abc'
            
        % artificial bee colony optimsation
        %------------------------------------------------------------------
        fprintf('Performing ABC Optimisation\n');     
        [Px,BestCost] = atcm.optim.abcAS(@optimi,Px,(~~Px)/8)
        
        case 'fminsearch'
            
        % fminsearch
        %------------------------------------------------------------------
        fprintf('Performing fminsearch Optimisation\n');
        [X,F(i)] = fminsearch(@optimi,Px);

        case 'sa'
            
        % simulated annealing
        %------------------------------------------------------------------
        fprintf('Performing Simulated-Annealing Optimisation\n');

        % generate bounds as percentage around init
        LB = (Px-1);
        UB = (Px+1);        
        [X,F(i)] = simulannealbnd(@optimi,Px,LB,UB);
        
        case 'ga'
            
        % genetic algorithm
        %------------------------------------------------------------------
        fprintf('Performing Genetic-Algorithm Optimisation\n');
        [X,F(i)] = ga(@optimi,length(Px));
        
        case 'hj'
            
        % Hooke-Jeeves Algorithm
        %------------------------------------------------------------------
        fprintf('Performing Hooke-Jeeves Optimisation\n');
        [X,F(i),Iters] = atcm.optim.hookejeeves(nc(i), Px, repmat(0.1,[nc,1]), repmat(0.005,[nc,1]), [], 100, @optimi);

        case{'dcm','free_energy','nlsi','fe'}
            
        % Free energy minisation - default DCM routine
        %------------------------------------------------------------------
        fprintf('Performing Free-energy Optimisation\n');
        DCM2       = DCM;
        DCM2.M.pE  = sparse(Px');
        DCM2.M.pC  = diag(sparse( (~~Px)*0.125'));
        DCM2.M.p0E = P;
        DCM2.M.IS  = @atcm.optim.aINTw_ica;
        DCM2.M.iW  = cm{i};
        DCM2.M.C   = cx{i};
        DCM2.M.f0  = DCM.M.IS; 
        [Qp,Cp,Eh,F(i)] = atcm.optim.spm_nlsi_GN_as(DCM2.M,DCM2.xU,DCM2.xY);
        X = Qp;
        
        case{'particle'}
            
        % Particle swarm optimisation
        %------------------------------------------------------------------
        fprintf('Performing Particle-Swarm Optimisation\n');
        %[X] = particleswarm(@optimi,length(Px));       
        options = optimoptions('particleswarm');
        options.InitialSwarmMatrix = Px';
        [X] = particleswarm(@optimi,length(Px),Px-5,Px+5,options);
        
        case{'pattern'}
            
        % Pattern search optimisation
        %------------------------------------------------------------------
        fprintf('Performing Pattern-Search Optimsation\n');
        [X,F(i)] = patternsearch(@optimi, (Px));
        
        case{'ls'}
            
        % Laplacian sampling
        %------------------------------------------------------------------
        fprintf('Performing Free-energy Optimisation Using Laplacian Sampling\n');
        DCM2       = DCM;
        DCM2.M.pE  = sparse(Px');
        DCM2.M.pC  = diag(sparse( (~~Px)*0.125'));
        DCM2.M.p0E = P;
        DCM2.M.IS  = @atcm.optim.aINTw_ica;
        DCM2.M.iW  = cm{i};
        DCM2.M.C   = cx{i};
        DCM2.M.f0  = DCM.M.IS; 
        DCM2.M.FS  = 'spm_fs_csd';
        
        [Ep,qC,qh,F] = atcm.optim.asample(DCM2.M,DCM2.xU,DCM2.xY);
        X = Ep;
        
                        
        case{'surrogate'}
            
        % surrogate model optimisation
        %------------------------------------------------------------------
        fprintf('Performing Surrogate-Model Optimisation\n');
        [X,F(i)] = surrogateopt(@optimi,Px*-5,Px*5);
        
        
    end
    
    close
        
    % Compute real-world parameter values & Update for repeat sequence
    %----------------------------------------------------------------------
    X0 = X'*cm{i};
    P  = spm_unvec(X0,P);
    
end


end

function E = optimi(x,varargin)
global DD

% The Objective function with mapping between reduced set and full parameter
% vector

% Precision
Q = 1;

% Use DCM-specified precision matrix if available
if isfield(DD.xY,'Q') && any(DD.xY.Q(:))
    Q = DD.xY.Q;
end

% Compute the current model estimate, given parameters
pE  = x;                                 % real world parameter vector
iW  = DD.cm;                             % mapping from subspace to real
PP  = spm_unvec( (pE*iW).*spm_vec(DD.M.pE) , DD.M.pE);

% Minimise free energy:
%--------------------------------------------------------------------------
DM           = DD;
DM.M.Nmax    = 1;
DM.M.nograph = 1;
DM.M.pC      = DD.M.pC;
%E            = atcm.optim.fe0(DM.M,DM.xU,DM.xY,PP);

% Minimise squared error term
%--------------------------------------------------------------------------
IS = spm_funcheck(DD.M.IS);
yy = IS(PP,DD.M,DD.xU);
Y  = (DD.xY.y); 
try
    E  = sum( Q*(spm_vec(Y) - spm_vec(yy)).^2 );
catch
    E  = sum( (spm_vec(Y) - spm_vec(yy)).^2 );
end

% Plot current prediction point (CSD)
%--------------------------------------------------------------------------
PCSD(Y{1},yy{1},[],DD.M.Hz);
drawnow;
hold off;


end


function [Px, cm, j] = clusteroptimse_maxi(P,DCM,nc)

% Parameters
V  = spm_vec(DCM.M.pC);
Ep = spm_vec(P);
ip = ~~V;

% Functions
IS = spm_funcheck(DCM.M.IS);

% Jacobian for active parameter space
j  = jaco(IS,Ep,P,DCM.M,DCM.xU,DCM.xY,ip);

% find parameters with biggest effects
for k = 1:nc
    Mj        = max(j');
    [~,p(k)]  = max(Mj);
    j(p(k),:) = 0;
end

cm    = zeros(nc, size(j,1));
for i = 1:nc
    cm( i, p(i) ) = 1;
end


% Make param vector for reduced space: [diag(Px)*cm]*P
Px = ones(size(cm,1),1);

end

function [Px, cm, cx, j] = clusteroptimse0(P,DCM,nc,CovP)

% Parameters
V  = spm_vec(DCM.M.pC);
Ep = spm_vec(P);
ip = ~~V;

% Functions
IS = spm_funcheck(DCM.M.IS);

if nc ~= length(find(ip));

    % Jacobian for active parameter space
    j  = jaco(IS,Ep,P,DCM.M,DCM.xU,DCM.xY,ip);
    j(~V,:) = NaN;
    
    % If parameter covariance supplied, incorporate
    try
        subj = CovP(~~V,~~V).*j(~~V,:);
        fprintf('Including prior estimation of param covariance in clustering\n');
    catch
        % only cluster jaco elements for v(P) > 0
        subj = j(~~V,:);
    end
    
    %j(isnan(j)) = 0; 
    
    % only cluster jaco elements for v(P) > 0
    subj = j(~~V,:);
    
    % Hierarchically cluster effects
    [C,iW] = clusters(subj,nc);
    
    % Put sub-param matrix into full size matrix
    ncm = zeros(size(iW,2),length(V));
    ncm(:,~~V) = iW'; 
    
    cm = ncm;
    cx = C;
    
    % Make param vector for reduced space: [diag(Px)*cm]*P
    Px = ones(size(cm,1),1);
else
    fprintf('No reduction selected\n');
    j  = 0;
    %cm = eye(nc);
    cm = zeros(length(find(ip)),length(ip));
    cm(:,find(ip)) = eye(nc);
    Px = ones(nc,1);
end

end

function [C,iW] = clusters(j,nc)

fprintf('Hierarchical Clustering\n');

% Hierarchical clustering on jaco covariance
%CV = cov(j');
CV = TSNorm(j,1,[],1);
%CV = j;
CV(isnan(CV)) = 0;

% Fixed dimensionality if num comps not specified
if isempty(nc)
    nc = round( rank(CV) /2);
end

% % The clustering routine
% c = clusterdata(real(CV),'linkage','ward','maxclust',nc);
% 
% % Generate cluster matrix - maps components to (real) parameters
% cm    = zeros(nc, size(j,1));
% for i = 1:nc
%     cm( i, find(c==i)) = 1;
% end


% OR use an SVD reduction
[C, A, W] = fastica(full(real(CV)),'numOfIC',nc);
nc = size(C,1);
iW = pinv(W);

% use: (C(n,:)'*iW(:,n)')';


end

function j = jaco(IS,Ep,P,M,U,Y,ip)

% Compute the Jacobian matrix using variable step-size
n  = 0;
%fx = spm_vec(IS(spm_unvec(Ep,P),M,U));
warning off ;

for i = 1:length(Ep)
    if ip(i)
        
        % Print progress
        n = n + 1;
        if n > 1; fprintf(repmat('\b',[1,length(str)])); end
        str  = sprintf('Computing Gradients [ip %d / %d]',n,length(find(ip)));
        fprintf(str);
        
        % Compute Jacobi: A(j,:) = ( f(x+h) - f(x-h) ) / (2 * h)
        P0     = Ep;
        P1     = Ep;
        d      = P0(i) * 0.01;
        P0(i)  = P0(i) + d  ;
        P1(i)  = P1(i) - d  ;
        
        %f0     = spm_vec(IS(spm_unvec(P0,P),M,U));
        %f1     = spm_vec(IS(spm_unvec(P1,P),M,U));
        
        % derivs with respect to error
        f0     = objectf(IS,spm_unvec(P0,P),M,U,Y);
        f1     = objectf(IS,spm_unvec(P1,P),M,U,Y);
        j(i,:) = (f0 - f1) / (2 * d);
        
        % Alternatively, include curvature
        %deriv1 = (f0 - f1) / 2 / d;
        %deriv2 = (f0 - 2 * fx + f1) / d ^ 2;
        %j(i,:) = deriv1 ./ deriv2;
    end
end

warning on;
fprintf('\n');
end

function e = objectf(IS,P,M,U,Y)

% An objective function
y = spm_vec( IS(P,M,U) );
X = spm_vec(Y.y);
e = X - y;

end

function X = makef(w,Fq,Amp,Wid)
%
%
%
% AS


try Wid ; catch Wid = 2; end
try Amp ; catch Amp = 2; end


X   = 0*w;
f   = findthenearest(Fq,w);
f   = f(1);
X   = X + Amp * exp( -(w-f).^2 / (2*(2*Wid)^2) );

end