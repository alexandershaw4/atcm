function [Q] = spm_gen_Q_as(P,X)
% Helper routine for spm_gen routines
% FORMAT [Q] = spm_gen_Q(P,X)
%
% P - parameters
% X - vector of between trial effects
% c - trial in question
%
% Q - trial or condition-specific parameters
%
% This routine computes the parameters of a DCM for a given trial, where
% trial-specific effects are deployed according to a design vector X. The
% parameterisation follows a standard naming protocol where, for example,
% X(1)*P.B{1} + X(2)*P.B{2}... adjusts P.A for all (input) effects encoded
% in P.B.
% P.BN and P.AN operate at NMDA receptors along extrinsic connections
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_gen_Q.m 5922 2014-03-18 20:10:17Z karl $


% condition or trial specific parameters
%==========================================================================
Q  = P;


% % trial-specific effects on C (first effect only)
% %--------------------------------------------------------------------------
% try
%     Q.C = Q.C(:,:,1) + X(1)*P.C(:,:,2);
% end

% trial-specific effects on A (connections)
%--------------------------------------------------------------------------
for i = 1:length(X)
    
    % Alex: allow full trial specific inputs [e.g. roving]
    % trial-specific effects on C (first effect only)
    %--------------------------------------------------------------------------
    try
        if i > 1; Q.C(:,:,i) = Q.C(:,:,i) + X(i)*P.C(:,:,i-1); end
    end

    
    % extrinsic (driving) connections
    %----------------------------------------------------------------------
    for j = 1:length(Q.A)
        
        Q.A{j} = Q.A{j} + X(i)*P.B{i};
        
        % CMM-NMDA specific modulation on extrinsic NMDA connections
        %------------------------------------------------------------------
        try
            Q.AN{j} = Q.AN{j} + X(i)*P.BN{i};
        end
        
    end
    
    % modulatory connections
    %----------------------------------------------------------------------
    try
        Q.M  = Q.M + X(i)*P.N{i};
    end
    
%     % intrinsic connections
%     %----------------------------------------------------------------------
%     if ndims(P.G)==2
%         try        
%             Q.G = P.G * X(i); 
%         end
%     elseif ndims(P.G) == 3
%         % For Rachael's LTP model, allow 3D G
%         % i.e. a full connectivity matrix for each trial
%         try
%             Q.G = P.G(:,:,i) * X(i);
%         end
%     end
%     
%     if isfield(P,'T1')
%         % trial specific modulation of TCs
%         Q.T1 = P.T1(i,:) * X(i);
%     end

    if isfield(P,'G')
        try
            Q.H      = Q.H + X(i) * P.G(:,:,i);
        end
    end
    if isfield(P,'T1')
        try
            Q.T(1:2) = Q.T(1:2) + X(i) * P.T1(i,:);
        end
    end
        
    
end