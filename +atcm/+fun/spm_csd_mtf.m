function [y,w,s,g] = spm_csd_mtf(P,M,U)
% Spectral response of a NMM (transfer function x noise spectrum)
% FORMAT [y,w,s,g] = spm_csd_mtf(P,M,U)
% P - parameters
% M - neural mass model structure
% U - trial-specific effects (induces expansion around steady state)
%
% y - {y(N,nc,nc}} - cross-spectral density for nc channels {trials}
%                  - for N frequencies in M.Hz [default 1:64Hz]
% w - frequencies
% s - modulation transfer functions (complex)
% g - normalised modulation transfer function (true Granger causality)
%
% When called with U this function will return a cross-spectral response
% for each of the condition-specific parameters specified in U.X; otherwise
% it returns the complex CSD for the parameters in P (using the expansion
% point supplied in M.x)
%
% When the observer function M.g is specified the CSD response is
% supplemented with channel noise in sensor space; otherwise the CSD
% pertains to hidden states.
%
% NB: requires M.u to specify the number of endogenous inputs
% This routine and will solve for the (hidden) steady state and use it as
% the expansion point for subsequent linear systems analysis (if trial
% specific effects are specified).
%
% See also:
%  spm_ccf2csd.m, spm_ccf2mar, spm_csd2ccf.m, spm_csd2mar.m, spm_mar2csd.m,
%  spm_csd2coh.m, spm_dcm_mtf.m, spm_Q.m, spm_mar.m and spm_mar_spectral.m
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
%
% Karl Friston
% $Id: spm_csd_mtf.m 6112 2014-07-21 09:39:53Z karl $
%
%
% hacked by as to incorporate into +atcm tools
%

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
if isfield(M,'Hz')
    w    = M.Hz;
else
    w    = 1:64;
    M.Hz = w;
end


% Alex additions

% w, initial states, dt | fs is specified & generate sampletimes & inputs
%--------------------------------------------------------------------------
w     = M.Hz;                     % FoI (w)
x     = M.x;                      % model (hidden) states
dt    = 1/1200;                   % hard-wired 400 hz
Fs    = 1/dt;                     % sampling frequency
tn    = 2;                        % sample window length, in seconds
pst   = 1000*((0:dt:tn-dt)');     % peristim times we'll sample

% For constant (DC) inputs...
%------------------------------------------------------------------
mu    = exp(P.R(1));              % mean amplitude
drive = ones(length(pst),1)*mu;   % amplitude (constant) over time
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
        %M.x   = atcm.fun.solvefixedpoint(Q,M);
    end
    
    % transfer functions (FFT of kernel)
    %----------------------------------------------------------------------
    S     = spm_dcm_mtf(Q,M);
    
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

                DoNoiseGamma = 1;
                if DoNoiseGamma
                    Gam = atcm.fun.makef(w,exp(P.G(1))*60,exp(P.G(2))*.05,4);
                    G(:,i,j) = G(:,i,j) + Gam(:);
                end
                DoSmoothIntro = 0;
                if DoSmoothIntro
                    G(1:3,i,j) = smooth( G(1:3,i,j) ); 
                end                
                
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



