function [X0,X,cm,F,j,others] = clusteroptimse(P,DCM,method,nc)
global DD
% Optimisation / inversion of DCM models using parameter space partitioning
% and optimisation via minimising either:
%   e = sum( data - model ).^2 
%
% or the free-energy:
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
% - artificial bee colony (ABC) - sampling based
% - fminsearch (gradient based)
% - hooke-jeeves pattern search algorithm (gradient free)
% - simulated annealing (gradient free, requires global optimisation toolbox)
% - genetic algorithm (requires global optimisation toolbox)
% - particle swarm
% - patternsearch
% - laplacian sampling (see spm_nlsi_LS)
% - surrogate model optimisation (requires global optimisation toolbox)
% - neural network algorithm
% - unscented Kalman filter
% - Coyote optimsation algorithm
% - Cuckoo optimsation algorithm
% - ADAM - gradient descent with adaptive learning rate
% - Bayesian optimsation
% - Bayesopt (Matlab)
% - Facebook's Nevergrad (Python) gradient free optimiser
% - direct: [di]viding [rect]angles
% - BADS - Bayesian adaptive direct search
% - VBMC - Variational Bayesian Monte Carlo
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
%     'NNA'       - neural network algorithm optimser
%     'ukf'       - uscented Kalman filter optimsation
%     'COA'       - Coyote optimsation algorithm
%     'cuckoo'    - Cuckoo optimsation algorithm
%     'ADAM'      - gradient descent with adaptive learning rate
%     'bo'        - Bayesian optimsation
%     'bayesopt'  - Bayesopt (Matlab)
%     'nevergrad' - Facebook's (Python) gradient free optimiser
%
% nc     = number of components (params) to reduce the full param set to.
% (nc can be a vector to iterate over)
%
% AS2018

DD   = DCM;
DD.P = P;
DD.Gmax = -inf;

others = [];

sf = 1; % scale factor on variances - i.e. bigger = wider var

% % force replace q
% q      = spm_Q(1/2,77,1); 
% q      = q*diag(4:80)*q;
% DCM.xY.Q = q;

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
    [Px, cm{i}, j{i}] = fs(P,DCM,nc(i),pC);
    DD.cm             = cm{i};
    DD.doplot         = 1;
    
    % Optimise the model.... options:
    %----------------------------------------------------------------------
    switch method
        
        case 'abc'
            
        % artificial bee colony optimsation
        %------------------------------------------------------------------
        fprintf('Performing ABC Optimisation\n');    
        c   = constr(Px');
        DD.doplot = 0;
        [X,F(i)]  = atcm.optim.abcAS(@optimi,Px,c*sf);
        
        case 'fminsearch'
            
        % fminsearch
        %------------------------------------------------------------------
        fprintf('Performing fminsearch Optimisation\n');
        [X,F(i)] = fminsearch(@optimi,Px);
        X = X(:);
        
        case 'sa'
            
        % simulated annealing
        %------------------------------------------------------------------
        fprintf('Performing Simulated-Annealing Optimisation\n');

        % generate bounds as percentage around init
        c   = constr(Px');
        LB  = (Px-(c'*sf));
        UB  = (Px+(c'*sf));  
        
        % options
        %saopts = optimoptions('simulannealbnd','MaxFunctionEvaluations', 1000);
        saopts = optimoptions('simulannealbnd','MaxTime', 86400/2); % 12hours
        
        [X,F(i),exitFlag,output] = simulannealbnd(@optimi,Px,LB,UB,saopts);
        others = output;
        
        case 'ga'
            
        % genetic algorithm
        %------------------------------------------------------------------
        fprintf('Performing Genetic-Algorithm Optimisation\n');
        c   = constr(Px') * 1e-2;
        LB  = (Px-(c'*sf));
        UB  = (Px+(c'*sf));  
        
        %[X,F(i)] = ga(@optimi,length(Px),[],[],[],[],LB,UB);
        gaopts = optimoptions('ga','MaxTime', 86400/2); % 12hours
        
        [X,F(i)] = ga(@optimi,length(Px),[],[],[],[],LB,UB,[],[],gaopts);        
        
        case 'hj'
            
        % Hooke-Jeeves Algorithm
        %------------------------------------------------------------------
        fprintf('Performing Hooke-Jeeves Optimisation\n');
        [X,F(i),Iters] = atcm.optim.hookejeeves(nc(i), Px, Px/10, Px-1, [], 100, @optimi);

        case{'dcm','free_energy','nlsi','fe'}
            
        % Free energy minisation - default DCM routine
        %------------------------------------------------------------------
        fprintf('Performing Free-energy Optimisation\n');
        DCM2       = DCM;
        DCM2.M.pE  = sparse(Px');
        DCM2.M.pC  = diag(sparse( (~~Px)*0.125'));
        DCM2.M.p0E = P;
        DCM2.M.IS  = @atcm.optim.aINTw_reduce;
        DCM2.M.cm  = cm{i};
        DCM2.M.f0  = DCM.M.IS; 
        [Qp,Cp,Eh,F(i)] = atcm.optim.spm_nlsi_GN_as(DCM2.M,DCM2.xU,DCM2.xY);
        X = Qp;
        
        case{'particle'}
            
        % Particle swarm optimisation
        %------------------------------------------------------------------
        fprintf('Performing Particle-Swarm Optimisation\n');
        %[X] = particleswarm(@optimi,length(Px));
        c   = constr(Px');
        LB  = (Px-(c'*sf));
        UB  = (Px+(c'*sf));  
        options = optimoptions('particleswarm','MaxTime', 86400/2);
        options.InitialSwarmMatrix = Px';
        [X,F(i)] = particleswarm(@optimi,length(Px),LB,UB,options);
        
        case{'pattern'}
            
        % Pattern search optimisation
        %------------------------------------------------------------------
        fprintf('Performing Pattern-Search Optimsation\n');
        c   = constr(Px');
        LB  = (Px-(c'*sf));
        UB  = (Px+(c'*sf));  
        [X,F(i)] = patternsearch(@optimi, Px,[],[],[],[],LB,UB);
        
        case{'ls'}
            
        % Laplacian sampling
        %------------------------------------------------------------------
        fprintf('Performing Free-energy Optimisation Using Laplacian Sampling\n');
        DCM2       = DCM;
        DCM2.M.pE  = sparse(Px');
        DCM2.M.pC  = diag(sparse( (~~Px)*0.125'));
        DCM2.M.p0E = P;
        DCM2.M.IS  = @atcm.optim.aINTw_reduce;
        DCM2.M.cm  = cm{i};
        DCM2.M.f0  = DCM.M.IS; 
        DCM2.M.FS  = 'spm_fs_csd';
        
        [Ep,qC,qh,F] = atcm.optim.asample(DCM2.M,DCM2.xU,DCM2.xY);
        X = Ep;
        
        
        case{'ls_long'}
            
        % Laplacian sampling
        %------------------------------------------------------------------
        fprintf('Performing Free-energy Optimisation using [Long] Laplacian Sampling\n');
        DCM2       = DCM;
        DCM2.M.pE  = sparse(Px');
        DCM2.M.pC  = diag(sparse( (~~Px)*0.125'));
        DCM2.M.p0E = P;
        DCM2.M.IS  = @atcm.optim.aINTw_reduce;
        DCM2.M.cm  = cm{i};
        DCM2.M.f0  = DCM.M.IS; 
        DCM2.M.FS  = 'spm_fs_csd';
        
        [Ep,qC,qh,F] = spm_nlsi_LS_long(DCM2.M,DCM2.xU,DCM2.xY);
        X = Ep;        
        
        case{'ls_gauss'}
            
        % Laplacian sampling
        %------------------------------------------------------------------
        fprintf('Performing Free-energy Optimisation using [Long] Laplacian Sampling\n');
        DCM2       = DCM;
        DCM2.M.pE  = sparse(Px');
        DCM2.M.pC  = diag(sparse( (~~Px)*0.125'));
        DCM2.M.p0E = P;
        DCM2.M.IS  = @atcm.optim.aINTw_reduce;
        DCM2.M.cm  = cm{i};
        DCM2.M.f0  = DCM.M.IS; 
        DCM2.M.FS  = 'spm_fs_csd';
        
        [Ep,qC,qh,F] = atcm.optim.spm_nlsi_LS_GausAS(DCM2.M,DCM2.xU,DCM2.xY);
        X = Ep;           
                
        case{'surrogate'}
            
        % surrogate model optimisation
        %------------------------------------------------------------------
        fprintf('Performing Surrogate-Model Optimisation\n');
        c   = constr(Px');
        LB  = (Px-(c'*sf));
        UB  = (Px+(c'*sf));  
        DD.doplot = 0;
        opts = optimoptions('surrogateopt','PlotFcn','surrogateoptplot');
        [X,F(i)] = surrogateopt(@optimi,LB,UB,opts);
        
        case 'NNA'
            
        % Neural network algorithm
        %------------------------------------------------------------------
        fprintf('Performing Neural-Network Algorithm Optimisation\n');
        c   = constr(Px');
        LB  = (Px-(c'*sf));
        UB  = (Px+(c'*sf));  

        [X,F(i)]=NNA_Const_Ver1(@optimi,@constr,LB',UB',length(Px));
        
        case{'ukf','ekf'};
        
        % unscented kalman filter optimser
        %------------------------------------------------------------------
        fprintf('Performing Unscented Kalman Filter Optimisation\n');
        c     = constr(Px');
        %[X,F] = ekfopt(@optimi,Px(:),1e-4,diag(c),eye(length(c)));
        [X,F] = ekfopt(@optimi,Px(:),1e-4);
        %[X,F] = ekfopt(@fakeDM,Px(:),1e-4);
        
        
        case 'COA'
            
        % Coyote Optimization Algorithm
        %------------------------------------------------------------------
        fprintf('Performing Coyote Optimization Algorithm\n');
        c           = constr(Px');
        dx          = ( Px.*(c') );
        D           = 77;                     % Problem dimension
        lu          = [(Px-(dx*sf))';(Px+(dx*sf))']; % Seach space
        nfevalMAX   = 1000*D;                % Stopping criteria
        Np          = 10;                     % Number of packs
        Nc          = 10;                     % Number of coyotes
        [X,F] = COA(@optimi, lu, nfevalMAX,Np,Nc);

        case 'cuckoo'
        
        % Cuckoo optimsation algorithm
        %------------------------------------------------------------------
        fprintf('Performing Cuckoo Optimization Algorithm\n');
        c   = constr(Px');
        LB  = (Px-(c'));
        UB  = (Px+(c'));  
        [X,F] = cuckooMain(@optimi,mean(LB),mean(UB),length(c));
        
        case 'adam'
        
        % ADAM optimsation algorithm
        %------------------------------------------------------------------
        % gradient descent with Adaptive learning rates
        fprintf('Performing GD w/ ADAM optimsation\n');
        [X, F] = fmin_adam(@optimi, Px');
        
        case 'bo'
            
        % Bayesian optimsation algorithm
        %------------------------------------------------------------------
        % IMGPO with a default setting used in a NIPS paper
        % http://lis.csail.mit.edu/code/imgpo.html
        % K. Kawaguchi, L. P. Kaelbling, T. Lozano-Pérez. 
        % Bayesian Optimization with Exponential Convergence.
        % In Advances in Neural Information Processing (NIPS), 2015.
        fprintf('Performing Bayesian Optimsation w/ Exp Convergence\n');
        c   = constr(Px');
        LB  = (Px-(c'*sf));
        UB  = (Px+(c'*sf));
        x_input_domain = [LB UB];
        
        result_diplay  = 1;
        result_save    = 1;
        plot_func      = 0;
        plot_point     = 1;
        DD.doplot      = 0;
        
        [X,F] = IMGPO_default_run(@optimi, 1, x_input_domain, 128, ...
                     result_diplay, result_save, plot_func, plot_point);
        
              
                 
        case 'bayesopt'
            
        % Bayesian optimsation algorithm
        %------------------------------------------------------------------
        fprintf('Performing (Matlab) Bayesian Optimsation\n');
        c   = constr(Px');
        LB  = (Px-(c'*sf));
        UB  = (Px+(c'*sf));
        
        for ip = 1:length(Px)
            name = sprintf('Par%d',ip);
            xvar(ip) = optimizableVariable(name,[LB(ip) UB(ip)],'Optimize',true);
        end

        reps    = 128;
        explore = 0.2;
        RESULTS = bayesopt(@optimi,xvar,'IsObjectiveDeterministic',true,...
                        'ExplorationRatio',explore,'MaxObjectiveEvaluations',reps,...
                        'AcquisitionFunctionName','expected-improvement-plus');

        % Best Actually observed model
        MinErr = RESULTS.MinObjective;
        F   = RESULTS.XAtMinObjective;
        X   = Fit.Variables;
        
        
        case 'direct'
        
        % DIRECT (DIviding RECTangles) optimsation algorithm
        %------------------------------------------------------------------
        fprintf('Performing DIRECT Optimization Algorithm\n');
        c   = constr(Px');
        LB  = (Px-(c'*sf));
        UB  = (Px+(c'*sf));
        bounds = [LB UB];

        %    We tell DIRECT that the globalmin = 3
        %    It will stop within 0.01% of solution
        options.testflag  = 1; 
        options.globalmin = 0; 
        options.showits   = 1;
        options.tol       = 0.01;
        options.maxits    = 10;
        
        % Pass Function as part of a Matlab Structure
        Problem.f = @optimi;
        
        % 3. Call DIRECT
        [F,X,hist] = Direct(Problem,bounds,options);

        
        
        case 'nevergrad'
        % Doesn't currently support assing matlab function handle :(
        %
        % Bayesian optimsation algorithm
        %------------------------------------------------------------------
        fprintf(['Facebook Nevergrad Optimsation doesnt yet support ',...
            'passing matlab function handles\n']);
        return;

        %pyversion /Users/Alex/miniconda3/envs/py6/bin/python
        opt = py.importlib.import_module('nevergrad.optimization');

        %optimizer = opt.optimizerlib.OnePlusOne(Px', 100);
        inst = py.importlib.import_module('nevergrad.instrumentation');

        %instrum = inst.Instrumentation(inst.var.Array(length(Px)), inst.var.Array(1).asfloat())
        instrum = inst.Instrumentation(inst.var.Array(1,length(Px)))
        optimizer = opt.optimizerlib.OnePlusOne(instrum, 100)

        recommendation = optimizer.optimize(@optimi)

        case {'aoptim','ao','aoptimgrad'}
            
        % Alex's gradient + sampling optimsation routine
        %------------------------------------------------------------------
        fprintf('Performing AO optimisation\n');
        DD.doplot = 0;
        c         = constr(Px');
        c         = c * sf;
        
        % - function is not an objective, but an integrator returning a vector
        % - supply vector target to aoptim, it cmputes its own object fun                
        %[X,F(i)]  = agradoptimi(@fakeDM,Px(:),c,DCM.xY.y{1});
        Q = DCM.xY.Q;
        Q = full(diag(Q));
        Q = smooth(Q);
        Q = diag(Q);
        
        %Q = diag( 1+(1:77)'/77 );;
        
        Q = diag(Q)*20;
        Q = diag(Q);
        
        %[X,F(i),Cp]  = AO(@fakeDM,Px(:),c,DCM.xY.y,256,12*4,[],1e-3);
        
        %[X,F(i),Cp]  = AO(@fakeDM,Px(:),c,DCM.xY.y,64*4,12*4,[],1e-3,0.000001);
          
        
        [X,F(i),Cp,Pp,History]  = AO(@fakeDM,Px(:),c,DCM.xY.y,32,12*4,[],...
                                    -inf,1e-12,2,0,'fe',0,1,[],1);
        
        
        %[X,F(i),Cp]  = AO(@optimi,Px(:),c,[],64*4,12*4,[],1e-3,0.000001);
        
        
        %[X,F(i),Cp]  = AO(@fakeDM,Px(:),c,DCM.xY.y,64,12*4,Q,1e-6,0.000001);

        % return non-reduced output - the unpacking is different for this
        % option:
        x0  = (diag(X) * cm{i});           
        x0  = sum(x0);
        x0(x0==0) = 1;
        x0(isnan(x0)) = 1;

        x1  = x0.*spm_vec(P)';        % P not DCM.M.pE!
        PP  = spm_unvec(x1,DCM.M.pE);        
        
        X0 = spm_vec(PP);
        return;
        
        
        case 'bads'
            
        % Bayesian adaptive direct search
        %------------------------------------------------------------------
        fprintf('Performing Bayesian adaptive direct search optimisation\n');
        c   = constr(Px');
        LB  = (Px-(c'*sf));
        UB  = (Px+(c'*sf));
        
        [X,F(i)] = bads(@optimi,Px(:)',LB',UB',LB',UB');
        
        case 'vbmc'
            
        % Variational Bayesian Monte Carlo
        %------------------------------------------------------------------
        fprintf('Performing Variational Bayesian Monte Carlo optimisation\n');
        c   = constr(Px');
        sf  = 1e-3;
        LB  = (Px-(c'*sf));
        UB  = (Px+(c'*sf));
        
        
        % wrapper the objective to make it negative! (log-likelihood maximi)
        fx0 = @(x) -optimi(x);
        
        [vp,elbo,elbo_sd] = vbmc(fx0,Px(:)',[],[],LB',UB');
            
        X0 = vp;
        X  = elbo;
        clear cm j
        cm = elbo_sd;
        j = 0;
        
        return;
        
        case 'gn'
            
        % A Gauss-Newton method
        %------------------------------------------------------------------
        fprintf('Performing Gauss-Newton optimisation\n');
        params{2} = DD.xY.y{1};
        [X,FLAG,nIter,X,F(i),ALPHA] = gaussnewton(@fakeDM,Px,params);
                    
        
        case 'gpso'
            
        % Jonathan Hadida's gaussian process optimser
        %------------------------------------------------------------------
        fprintf('Performing GPMO optimisation\n');
        objfun = @optimi; % 
        
        c   = constr(Px');
        LB  = (Px-(c'*sf));
        UB  = (Px+(c'*sf));

        domain = [LB UB]; % domain of different size in each dimension

        obj = GPSO(); % create instance for optimisation
        Neval = 256; % budget of 50 evaluations
        output = obj.run( objfun, domain, Neval, 'tree', 3 ); % explore 
        % children leaf by growing partition tree 3 levels deep
        
        case 'none'
            
            % just return the jacobian matrix
            F  = 0;
            X  = Px; 
                
        case {'PR','minimise','hinton'};
            
            % The Polack-Ribiere flavour of conjugate gradient optimiser
            
            [X, F, ix] = PR_minimize(Px, @optimi, 0);
            
        case 'qpso'
            
            % quantum particle swarm opt
            %------------------------------------------------------------------
            fprintf('Performing quantum particle swarm optimisation\n');
            objfun = @optimi; % 

%             c   = constr(Px');
%             LB  = (Px-(c'*4));
%             UB  = (Px+(c'*4));
            
        case 'WOA'
            
            % whale optimsation algorithm
            %------------------------------------------------------------------
            fprintf('Performing Whale Optimisation Algorithm\n');
            
            SearchAgents_no=30; % Number of search agents
            Max_iteration=500; % Maximum numbef of iterations
            
            c   = constr(Px');
            LB  = (Px-(c'*sf));
            UB  = (Px+(c'*sf));
            
            [F,X,WOA_cg_curve]=WOA(SearchAgents_no,Max_iteration,LB',UB',length(LB),@optimi);

            
            
    end
    
    
    %close
    
    % Compute real-world parameter values
    %----------------------------------------------------------------------
    X0 = sum(diag(X)*cm{i});
    X0(X0==0) = 1;
    X0 = full(X0.*exp(spm_vec(P)'));
    X0 = log(X0);
    
    % Update for repeat sequence
    %----------------------------------------------------------------------
    P  = spm_unvec(X0,P);
end


end

% function [x,j] = AddJac(fun,x0)
% 
% x = fun(x0);
% V = ones(size(x0));
% [j,ip] = jaco(fun,x0,V);
% 
% 
% end


function [y,J] = fakeDM(Px,varargin)
global DD

pE  = DD.P;                          % real world parameter vector
hpE = DD.cm;                       % mapping from subspace to real
x0  = (diag(Px) * hpE);           % indexing
x0  = sum(x0);
x0(x0==0) = 1;
x0(isnan(x0)) = 1;

x1  = x0.*spm_vec(pE)';     % new (full) parameter values
PP  = spm_unvec(x1,pE);        % proper structure


M    = DD.M;
M.pC = ( ~~real(spm_vec(PP)) +0.125 );
M.pC = spm_unvec(M.pC,PP);

IS = spm_funcheck(DD.M.IS);       % Integrator
y  = IS(PP,M,DD.xU);               % Prediction
y  = spm_vec(y);

% convert to residuals
%y = y - spm_vec(DD.xY.y);
y = real(y);

% plot
if DD.doplot
    plot(DD.xY.Hz,DD.xY.y{1},':',DD.xY.Hz,y); drawnow;
end

if nargout > 1
        
    % compute jacobi
    %V = ones(size(Px))/1e-6;
    V = DD.cm*spm_vec(DD.M.pC);
    %[J,ip] = jaco1(@fakeDM,Px, V.*(hpE*spm_vec(pE)) );
    [J,ip] = jaco1(@fakeDM,Px, V  );
    
    J = J';
    return;
end


end



function c = constr(x)
global DD

pC = spm_vec(DD.M.pC);
cm = DD.cm;
c  = (pC'*cm').*x;

%c = x;
end

function [E,J] = optimi(x,varargin)
global DD

% The Objective function with mapping between reduced set and full parameter
% vector

Gmax = DD.Gmax;

% Bayesopt returns x as a table rather than vecotr: reshape!
if istable(x)
    x = x.Variables;
end

% Precision
Q = 1;

% Use DCM-specified precision matrix if available
if isfield(DD.xY,'Q') && any(DD.xY.Q(:))
    Q = DD.xY.Q;
end

% Compute the current model estimate, given parameters
pE  = DD.P;                              % real world parameter vector
hpE = DD.cm;                             % mapping from subspace to real
x0  = (diag(x) * hpE);                   % indexing
x0  = sum(x0);                           % actual values (as vector)
x0(x0==0) = 1;                           % 0==1 so as not to change fixed params
try;  x1  = x0.*exp(spm_vec(pE))';       % new (full) parameter values
catch;x1  = x0.*exp(spm_vec(pE)); 
end
PP  = spm_unvec(log(x1),pE);             % proper (full model) structure


% Minimise free energy:
%--------------------------------------------------------------------------
DM           = DD;
DM.M.Nmax    = 1;
DM.M.nograph = 1;
DM.M.pC      = DD.M.pC;
% E            = atcm.optim.fe0(DM.M,DM.xU,DM.xY,PP,Gmax);
% E            = -real(E);
% 
% if E < Gmax
%     Gmax = E;
%     DD.Gmax = Gmax;
% end

% % Minimise squared error term
% %--------------------------------------------------------------------------
IS = spm_funcheck(DD.M.IS);
try
    try 
        [yy,w,s,g,t,pst] = IS(PP,DD.M,DD.xU);
    catch
        [yy,w] = IS(PP,DD.M,DD.xU);
    end
catch
    yy = spm_unvec( (spm_vec(DD.xY.y)*0)+inf , DD.xY.y);
    w  = DD.xY.Hz;
end

Y  = (DD.xY.y); 
try
    E  = sum( Q*(spm_vec(Y) - real(spm_vec(yy))).^2 );
catch
    E  = sum( (spm_vec(Y) - real(spm_vec(yy))).^2 );
end

% Plot current prediction point (CSD)
%--------------------------------------------------------------------------
%PCSD(Y{1},yy{1},[],DD.M.Hz);
%drawnow;
%hold off;
if DD.doplot
    warning off
    
    if ndims(Y{1}) < 3
        try
            subplot(211); plot(w,Y{1},w,yy{1});
        end
        try
            subplot(212); plot(pst,squeeze(s{1}(1,:,1,:)));
        end
            drawnow;
            hold off;
    
    else
    
        if ndims(Y{1})==3
            PCSD(Y{1},yy{1},[],w);drawnow;
        end
    end
    
    warning on
end

if nargout > 1
    % compute jacobi
    V = ones(size(x));
    [J,ip] = jaco1(@optimi,x,V);
end

end

function [j,ip] = jaco1(fun,x0,V)
% Compute the Jacobian matrix (parameter gradients) for a model using:
%
% j(ip,:) = ( f(x(ip)+h)  - f(x(ip)-h) )  / (2 * h)
%
%
% AS2019
global DD

IS = fun;
P  = x0(:);

% if nargin == 3; ip = find(V(:));
% else;           ip = 1:length(x0);
% end

if nargin == 3; ip = (V(:));
else;           ip = 1:length(x0);
end

j  = jacf(IS,P,ip);

j(isnan(j)) = 0;

end



function j = jacf(IS,P,ip)
global DD

% Compute the Jacobian matrix using variable step-size
n  = 0;
warning off ;

f0    = feval(IS,P);
fx = f0;
j     = zeros(length(P),length(f0(:))); % n param x n output
for i = 1:length(P)
    if ~~ip(i)
        
        % Print progress
        n = n + 1;
        if n > 1; fprintf(repmat('\b',[1,length(str)])); end
        str  = sprintf('Computing Gradients [ip %d / %d]',n,length(find(ip)));
        fprintf(str);
        
        % Compute Jacobi: A(j,:) = ( f(x+h) - f(x-h) ) / (2 * h)
        P0     = P;
        P1     = P;
        d      = P0(i) * ip(i);
        
        if d == 0
            d = 0.01;
        end
        
        P0(i)  = P0(i) + d  ;
        P1(i)  = P1(i) - d  ;
                
        f0     = (feval(IS,P0));
        f1     = (feval(IS,P1));
        j(i,:) = (f0 - f1) / (2 * d);
        
        % Alternatively, include curvature
        deriv1 = (f0 - f1) / 2 / d;
        deriv2 = (f0 - 2 * fx + f1) / d ^ 2;
        j(i,:) = deriv1 ./ deriv2;
    end
end

warning on;
fprintf('\n');
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

function [Px, cm, j] = clusteroptimse0(P,DCM,nc,CovP)

% Parameters
V  = spm_vec(DCM.M.pC);
Ep = spm_vec(P);
ip = ~~V;

% Functions
IS = spm_funcheck(DCM.M.IS);

if nc ~= length(find(ip));

    % Jacobian for active parameter space
    j  = jaco0(IS,Ep,P,DCM.M,DCM.xU,DCM.xY,ip);
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
    cm = clusters(subj,nc);
    
    % Put sub-param matrix into full size matrix
    ncm = zeros(nc,length(V));
    ncm(:,~~V) = cm; 
    
    cm = ncm;
    
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

function cm = clusters(j,nc)

fprintf('Hierarchical Clustering\n');

% Hierarchical clustering on jaco covariance
%CV = cov(j');
 CV = TSNorm(j,1,[],1);
 CV(isnan(CV)) = 0;

% CV = cdist(j,j);
% CV(isnan(CV)) = 0;

% Fixed dimensionality if num comps not specified
if isempty(nc)
    nc = round( rank(CV) /2);
end

% The clustering routine
c = clusterdata(real(CV),'linkage','ward','maxclust',nc);

% Generate cluster matrix - maps components to (real) parameters
cm    = zeros(nc, size(j,1));
for i = 1:nc
    cm( i, find(c==i)) = 1;
end

end

function j = jaco0(IS,Ep,P,M,U,Y,ip)

% Compute the Jacobian matrix using variable step-size
n  = 0;
fx = spm_vec(IS(spm_unvec(Ep,P),M,U));
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
        
        % derivs w.r.t error
        f0     = objectf(IS,spm_unvec(P0,P),M,U,Y);
        f1     = objectf(IS,spm_unvec(P1,P),M,U,Y);
        j(i,:) = (f0 - f1) / (2 * d);
        
        % Alternatively, include curvature
        deriv1 = (f0 - f1) / 2 / d;
        deriv2 = (f0 - 2 * fx + f1) / d ^ 2;
        j(i,:) = deriv1 ./ deriv2;
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