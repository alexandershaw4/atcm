function [X,F] = agradoptim(DCM,x0,V)
global aopt
% partially-gradient dependent optimisation routine, optimsised for DCMs
% but can be used for anything
%
% steps:
% (1)  compute parameter jacobian (dfdx) w.r.t error (see atcm.fun.jaco.m)
% (2)  feature selection: identify regions of biggest difference (precision weighted)
% (3)  use jaco to select candidate parameters to minimise error in identified region
% (4)  loop through identified parameters, evaluating sampling around the mean
% (5)  accept or reject samples
% (6)  return to step (2) onwards ... 
% (~7) occasionally recompute the jacobian
%
% If n-iterations pass without improvement, increase number of candidate
% parameters selected in (3)
%
%
% AS2019

% check functions
%--------------------------------------------------------------------------
aopt.DCM  = DCM;
x0        = full(spm_vec(x0));
V         = full(spm_vec(V));
np        = length(x0);
npi       = length(find(V));
w         = DCM.xY.Hz;
[e0,er,Q] = obj(x0);

% compute jacobian
%--------------------------------------------------------------------------
D0       = DCM;
D0.M.pC  = spm_unvec(V,DCM.M.pC);
[j,ip]   = atcm.fun.jaco(D0,spm_unvec(x0,DCM.M.pE));

% jac w.r.t. error
%--------------------------------------------------------------------------
je = j - repmat(er',[np,1]);

% convergence / tol & variables
%--------------------------------------------------------------------------
fprintf('Starting main optimsation routine...\n')
contol   = 1e-4; % convergence error point
nit      = 0;    % just a counter
innit    = 0;    % inner loop counter
MaxIt    = 256;  % max (outer) iterations
DidGood  = 0;    % flat for accepted updated
PauPnt   = round(linspace(1,MaxIt-25,25)); % pause points
NpC      = 8 ; % number of parameter to consider at each point
NpS      = 2 ; % number of parameter to find for each minimisation point
PlotGood = 1;  % flag to plot on acceppt on ondating
pmat     = ones(1,1,1,1)*inf; % parameter trajectory matrix (w/parm values)
emat     = ones(1,1,1,1)*inf; % parameter trajectory matrix (w/errors)
nac      = 0;  % counter for accepted updates

% initialise loop
while e0 > contol
    
    % update loop counters
    %----------------------------------------------------------------------
    nit          = nit + 1;
    innit        = 0;
    DidGood(nit) = 0;
    eAll(nit)    = e0;
    
    % compute regions of biggest error
    %----------------------------------------------------------------------
    ind = featsel(er,Q,NpC);
    
    % best n parameters for minimising this error
    %----------------------------------------------------------------------
    P0 = [];
    for nps     = 1:length(ind)
        [~,P0i] = atcm.fun.maxpoints(je(:,ind(nps)),NpS,'min');
        P0      = [P0(:);P0i(:)];
    end;P0      = unique(P0);
    
    % sample around the mean (+/-) for these parameters
    %----------------------------------------------------------------------
    for nip = 1:length(P0)
        px  = x0(P0(nip));
        vx  = V(P0(nip));
        
        % sampling: inner loop
        for k = 1:4
            innit = innit + 1;
            
            % sample from normal distribution
            pxk = px + vx*randn(1);
            dx  = x0;
            dx(P0(nip)) = pxk;
            
            % compute objective
            [de,der] = obj(dx);
           
            if de < e0
                e0 = de;  % new error
                er = der; % new difference
                x0 = dx;  % new parameters
                
                % print updates to screen
                pupdate(nit,innit,de,e0,'accept');
                DidGood(nit) = 1;
                nac          = nac + 1;
                
                % plot updated model prediction and errors
                if PlotGood
                   DoUpdatePlot(x0,e0,eAll); 
                end
                                
                continue;
            else
                % print & just try again
                pupdate(nit,innit,de,e0,'reject');
            end
    
        end    
    end
    
    % pause to recompute jacobian ocassionally
    %----------------------------------------------------------------------
    if ismember(nit,PauPnt)
        if nit > 10 && sum( DidGood(nit-10:nit)  ) > 0
            % recompute jaco
            [j,ip]   = atcm.fun.jaco(D0,spm_unvec(x0,DCM.M.pE));
            % jac w.r.t. error
            je = j - repmat(er',[np,1]);
        end
    end
    
    % if too long without imporovement, expand n params
    %----------------------------------------------------------------------
    if nit > 20 && sum( DidGood(nit-10:nit)  ) == 0
        fprintf('Increasing num parameters considered\n'); 
        NpC = NpC + 2;
    end
        
    % stop at max iteration
    %----------------------------------------------------------------------
    if nit == MaxIt
        X = x0;
        F = e0;
        return;
    end
    
    % convergence criterion
    %----------------------------------------------------------------------
    if e0 <= contol
        X = x0;
        F = e0;
        fprintf('Converged!\n');
        return;
    end
    
end

% If we don't converge, still return bests
X = x0;
F = e0;

end

function pupdate(it,innit,err,best,action)


fprintf('| Main It: %04i | Inner It: %04i | Err: %04i | Best: %04i | %s |\n',it,innit,err,best,action);


end



function DoUpdatePlot(x0,enow,eall)
global aopt

DD.DCM = aopt.DCM;

% Precision
%--------------------------------------------------------------------------
Q = 1;

% Use DCM-specified precision matrix if available
if isfield(DD.DCM.xY,'Q') && any(DD.DCM.xY.Q(:))
    Q = DD.DCM.xY.Q;
end

% Minimise squared error term
%--------------------------------------------------------------------------
IS = spm_funcheck(DD.DCM.M.IS);
P  = spm_unvec(x0,DD.DCM.M.pE);
FullPlot = 0;

try
    try
        [yy,w,s,g,t,pst] = IS(P,DD.DCM.M,DD.DCM.xU);
        FullPlot = 1;
    catch
        [yy,w] = IS(P,DD.DCM.M,DD.DCM.xU);
        FullPlot = 0;
    end
catch
    yy = spm_unvec( (spm_vec(DD.DCM.xY.y)*0)+inf , DD.DCM.xY.y);
    w  = DD.xY.Hz;
end

if FullPlot
    subplot(211); plot(w,DD.DCM.xY.y{1},':',w,yy{1},'r','linewidth',2);
    xlabel('Frequency (Hz)'); ylabel('PSD');title('Current Best Estimate');

    subplot(212); plot(1:length(eall),eall,'.','MarkerSize',20);hold on;
    plot(1:length(eall),eall,'k--');
    plot(length(eall)+1,enow,'.','MarkerSize',20); hold off;
    ylabel('Error.^2');xlabel('Improvments in error');title('Error');
    drawnow;
else
    plot(w,DD.DCM.xY.y{1},':',w,yy{1},'r','linewidth',2);
    xlabel('Frequency (Hz)'); ylabel('PSD');title('Current Best Estimate');
    drawnow;
end


end

    
function ind = featsel(er,Q,NpC)
% identify regions of w where the error in y is biggest

[~,ind] = atcm.fun.maxpoints(abs(Q*er),NpC);


end



function [e,er,Q] = obj(x0)
global aopt

DD.DCM = aopt.DCM;

% Precision
%--------------------------------------------------------------------------
Q = 1;

% Use DCM-specified precision matrix if available
if isfield(DD.DCM.xY,'Q') && any(DD.DCM.xY.Q(:))
    Q = DD.DCM.xY.Q;
end

% Minimise squared error term
%--------------------------------------------------------------------------
IS = spm_funcheck(DD.DCM.M.IS);
P  = spm_unvec(x0,DD.DCM.M.pE);

try
    try
        [yy,w,s,g,t,pst] = IS(P,DD.DCM.M,DD.DCM.xU);
    catch
        [yy,w] = IS(P,DD.DCM.M,DD.DCM.xU);
    end
    
catch
    yy = spm_unvec( (spm_vec(DD.DCM.xY.y)*0)+inf , DD.DCM.xY.y);
    w  = DD.xY.Hz;
end

% Error term
%--------------------------------------------------------------------------
Y  = (DD.DCM.xY.y); 
try
    e  = sum( Q*(spm_vec(Y) - spm_vec(yy)).^2 );
catch
    e  = sum( (spm_vec(Y) - spm_vec(yy)).^2 );
end

er = spm_vec(yy) - spm_vec(Y);

end

% function [pmat,emat] = updatelog(pmat,emat,ind,P0,pxk,e0,nit)
% % keep track of the parameter updates that minimise error - along with the
% % frequency they were selected to minimise w(ind), the value accepts (pxk)
% % and the global error (e0) so we can compute parameter trajectories
% 
% pmat(ind,P0,nit) = pxk;
% emat(ind,P0,nit) = e0;
% 
% end
% 
% function t = computetrajectories(pmat,emat)
% 
% nf  = size(pmat,1);
% np  = size(pmat,2);
% nit = size(pmat,3);
% 
% for i = 1:np
%     pm = squeeze(pmat(:,i,:));
%     em = squeeze(emat(:,i,:));
% 
%     xy   = [pm(:) em(:)];
%     p(i,:) = polyfit(xy(:,1),xy(:,2),1); % polyval(p{i},n)
%     
%     if any(isnan(p(i,:)))
%         p(i,:) = [0 0];
%     end
%     
% end
% 
% t = p;
% 
% end
