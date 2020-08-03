function DCM = parameters(DCM,Ns)
% collect initial (prior) parameters
%
% Need to rewrite this function by writing out the priors saved in matfile

%load('april23_priors.mat','pE','pC');

% load('april25_priors.mat','pE','pC');
% 
% pC.G = pC.G*0;
% 
% load('SomePrior','Pp');

if Ns == 1
    %load('May7.mat');

    %load NewDelayPriors.mat
    pth = fileparts(mfilename('fullpath'));
    %load([pth '/AugSpectralPriors']);
    load([pth '/NewModPriors']);
    
    DCM.M.pE = pE;
    DCM.M.pC = pC;
            
else
    
    % I need to put the priors from .mat above into here for multinode
    % setup - CBA right now
    fprintf('Priors for multi-node models are out of date\n');
    
    % NEEDS UPDATING FOR MULTINODE/REGION
    pth = fileparts(mfilename('fullpath'));
    pr = load([pth '/NewModPriors']);;
    
    % Extrinsics: restructure adjacency matrices
    %----------------------------------------------------------------------
    A     = DCM.A;
    A{1}  = A{1} | A{3};                              % forward
    A{2}  = A{2} | A{3};                              % backward

    % [F] SP -> SS & DP
    pE.A{1} = A{1}*32 - 32;
    pC.A{1} = A{1}/8;

    % [B] DP -> SP & SI
    pE.A{2} = A{2}*32 - 32; 
    pC.A{2} = A{2}/8;

    % [Extrinsic T->C] (none)
    pE.A{3} = A{1}*32 - 32; 
    pC.A{3} = A{1}/8;

    % [Extrinsic T->T] rt -> rt & rl
    AA      = A{1} | A{2};
    
    pE.A{4} = AA*32 - 32; 
    pC.A{4} = AA/8;
    
    % [Extrinsic T->T] rl -> rt & rl
    pE.A{5} = AA*32 - 32; 
    pC.A{5} = AA/8;

    % NMDA = same as AMPA
    for i = 1:length(pE.A)
        pE.AN{i} = pE.A{i};
        pC.AN{i} = pC.A{i};
    end

    % Beta's: input-dependent scaling
    %----------------------------------------------------------------------
    B     = DCM.B;
    for i = 1:length(B)
        B{i} = ~~B{i};
        pE.B{i} = B{i} - B{i};
        pC.B{i} = B{i}/8;
    end

    for i = 1:length(B)
        B{i} = ~~B{i};
        pE.BN{i} = B{i} - B{i};
        pC.BN{i} = B{i}/8;
    end
    
    % Intrinsic (local) parameters, per region
    %----------------------------------------------------------------------
    ns = length(A{1}); % number of regions / nodes
    np = 8;            % number of populations per region
    nk = 7;            % number of states per population
    
    % exogenous inputs
    %----------------------------------------------------------------------
    C     = DCM.C;
    C     = ~~C;
    pE.C  = C*32 - 32;
    pC.C  = C/8;
    
    % new multi-input model
    %pE.C = repmat(pE.C,[1 4]);
    %pC.C = repmat(pC.C,[1 4]);
    %pE.C = repmat([-0.1164 -0.0025 -0.0017 -0.0011],[ns 1]);
    
    %pE.C = repmat([-0.1641 -0.0023 -0.0013 -0.0028],[ns 1]);
    %pC.C = ~~pE.C/8;
    
    % Average Firing
    pE.S = pr.pE.S;
    %pE.S = [-0.1271 0.1873 -0.0859 -0.0357 -0.0286 0.2409 0.1187 -0.0886];
    pC.S = ones(1,8)*0.0625;
    
    % Average (Self) Excitation Per Population
    pE.G = zeros(ns,np);
    pC.G = zeros(ns,np); 
    %pC.G(1:2) = 1;
    
    % Average Membrane Capacitance
    %pE.CV = [-0.1148 0.0716 -0.0950 -0.1269 0.2081 0.1757 0.0406 -0.1243];
    pE.CV = pr.pE.CV;
    pC.CV = ones(1,8) * 0.0625;
    
    % Average Background Activity
    %pE.E = 0;
    pE.E = pr.pE.E;
    pC.E = 0;
    
    % Average Intrinsic Connectivity
    pE.H = repmat(pr.pE.H, [1 1 ns]);
    
    %pC.H = ~~pE.H / 8;
    pC.H =repmat([...
    0.0625         0    0.0625         0         0    0.0625         0    0.0625
    0.0625    0.0625    0.0625         0         0         0         0         0
    0.0625    0.0625    0.0625         0         0         0         0         0
         0    0.0625         0    0.0625    0.0625         0         0         0
         0         0         0    0.0625    0.0625         0         0         0
         0         0         0    0.0625    0.0625    0.0625         0    0.0625
         0         0         0         0         0         0    0.0625    0.0625
         0         0         0         0         0    0.0625    0.0625    0.0625], [1 1 ns]);
    
    % Receptor Time Constants (Channel Open Times)
    pE.T = repmat(pr.pE.T,[ns,1]);
    pC.T =  ( ones(ns,4) / 8);
    
    % Parameters on input bump: delay, scale and width
    %pE.R = [-0.0360 0 0];
    pE.R = pr.pE.R;
    pC.R = [0.25 0 0];
    
    % Delays: states - intrinsic & extrinsic
    %pE.D = [-0.0330 0.0034];
    
    pE.D = pr.pE.D;
    pC.D = [0 0]+1/8;
    
    
    
    %pE.D0 = [-0.0872 0.1470];
    pE.D0 = pr.pE.D0;
    pC.D0 = [1 1]/8;
    
    % Dipole position (not in use)
    pE.Lpos = sparse(3,0);
    pC.Lpos = sparse(3,0);
    
    % Electrode gain
    %pE.L = repmat(14.7361         ,[ns,1]);
    pE.L = repmat(pr.pE.L(1)         ,[ns,1]);
    pC.L = repmat(64              ,[ns,1]);
        
    % Contributing states
    %J           = zeros(1,np,nk) - 1000;
    %J(:,[1 2 4 6],1)    = log([.2 .8 .2 .2]);
    %J(isinf(J)) = -1000;
    %pE.J        = spm_vec(J)';
    %pC.J        = spm_vec(J)' * 0;
    
    %pE.J = [-0.2452 -0.2839 -0.2535 -0.2524];
    pE.J = pr.pE.J;
    pC.J = ones(1,4)/8;
        
    % Noise Components - a, b & c
    pE.a = repmat([-1.7311;0],[1 ns]);
    %pE.a = pr.pE.a;
    pC.a = repmat([1/16;0]   ,[1 ns]);
    
    pE.b = [-0.0664;0];
    %pE.b = pr.pE.b;
    pC.b = [1/16   ;0];
    
    pE.c = repmat([-24.4382;-0.3475],[1 ns]);
    pC.c = repmat([0;0],                      [1 ns]);
    
    % Neuronal innovations
%     pE.d = [    -1.0300
%    -0.2037
%     0.8057
%    -0.2735
%    -0.0366
%    -0.1578
%    -0.4520
%    -0.0749];

          
    pE.d = pr.pE.d;
    pC.d = ones(length(pE.d),1) / 8;
       
    pE.h = repmat(0,  [ns 1]);
    pC.h = repmat(0,  [ns 1]);
    
    pE.m = repmat(0,  [ns 1]);
    pC.m = repmat(0,  [ns,1]);
    
    pE.TV = [-0.0972 -0.7756 0.5177 0.0339 -0.0856 0.1213 -0.2327 -0.0248];
    pC.TV = ones(1,8)/16;
    
    pE.Mh = zeros(8,1);
    pC.Mh = zeros(8,1);
    
    pE.Hh = [0 0];
    pC.Hh = [0 0];
    
    pE.ID = [0.0301 0.0184 0.0120 -0.0286 0.0289 -0.0581 -0.0207 -0.0182];
    pC.ID = ones(1,8)/8;
    
    pE.Hn=pE.H;
    pC.Hn=pC.H;
    
    % Pack pE & pC into M
    DCM.M.pE = pE;
    DCM.M.pC = pC;
    
end
    
    
end
