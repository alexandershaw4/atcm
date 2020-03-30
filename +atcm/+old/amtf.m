function [y,w,s,g] = amtf(P,M,U)




% between-trial (experimental) inputs
%==========================================================================
try
    X = U.X;
    if ~size(X,1)
        X = sparse(1,0);
    end
catch
    
    % default inputs - one trial (no trial-specific effects)
    %----------------------------------------------------------------------
    X = sparse(1,0);
    
end


% compute log-spectral density
%==========================================================================

% frequencies of interest
%--------------------------------------------------------------------------
w     = M.Hz;
dt    = 1/1200;                   % hard-wired 400 hz
Fs    = 1/dt;                     % sampling frequency
tn    = 2;                        % sample window length, in seconds
pst   = 1000*((0:dt:tn-dt)');     % peristim times we'll sample
M.pst = pst;

% number of channels and exogenous (neuronal) inputs or sources
%--------------------------------------------------------------------------
nc   = M.l;
nw   = length(M.Hz);

% spectrum of innovations (Gu) and noise (Gs and Gn)
%--------------------------------------------------------------------------
if isfield(M,'g')
    [Gu,Gs,Gn] = spm_csd_mtf_gu(P,M.Hz);
else
    Gu         = spm_csd_mtf_gu(P,M.Hz);
    nc         = size(Gu,2);
end

% cycle over trials (experimental conditions)
%==========================================================================
for  c = 1:size(X,1)
    

    % condition-specific parameters
    %----------------------------------------------------------------------
    Q   = spm_gen_Q(P,X(c,:));
    
    % solve for steady-state - if exogenous inputs are specified
    %----------------------------------------------------------------------
    if nargin > 2
        M.x = spm_dcm_neural_x(Q,M);
    end
    
    % transfer functions (FFT of kernel)
    %----------------------------------------------------------------------
    S     = atcm.adcmmtf(Q,M);
    
    % predicted cross-spectral density
    %----------------------------------------------------------------------
    G     = zeros(nw,nc,nc);
    for i = 1:nw
        G(i,:,:) = sq(S(i,:,:))*diag(Gu(i,:))*sq(S(i,:,:))';
    end
    
    % save trial-specific frequencies of interest
    %----------------------------------------------------------------------
    g{c}  = G;
    s{c}  = S;
    
    
end

% and add channel noise
%==========================================================================
if isfield(M,'g')
    
    for c = 1:length(g)
        G = g{c};
        for i = 1:nc
            
            % channel specific noise
            %--------------------------------------------------------------
            try
                G(:,i,i) = G(:,i,i) + Gs(:,i);
            catch
                G(:,i,i) = G(:,i,i) + Gs(:,1);
            end
            
            % and cross-spectral density from common channel noise
            %--------------------------------------------------------------
            for j = 1:nc
                G(:,i,j) = G(:,i,j) + Gn;
            end
        end
        y{c} = G;
        
    end
else
    y = g;
end

% Granger causality (normalised transfer functions) if requested
%==========================================================================
if nargout > 3
    for c = 1:length(s)
        g{c} = spm_dtf2gew(s{c},Gu);
    end
end

% squeeze but ensure second dimension is returned as a common vector
%--------------------------------------------------------------------------
function [x] = sq(x)
if size(x,3) > 1, x = squeeze(x); else, x = x(:); end



