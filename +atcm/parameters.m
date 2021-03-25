function DCM = parameters(DCM,Ns,priorfile)
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

    if nargin < 3 || isempty(priorfile)
        %load NewDelayPriors.mat
        pth = fileparts(mfilename('fullpath'));
        %load([pth '/AugSpectralPriors']);
        load([pth '/+fun/Priors2021a']);
    else
        load(priorfile)
    end
    
    DCM.M.pE = pE;
    DCM.M.pC = pC;
            
else
    
    % I need to put the priors from .mat above into here for multinode
    % setup - CBA right now
    fprintf('Priors for multi-node models may be out of date\n');
    
        if nargin < 3 || isempty(priorfile)

        % NEEDS UPDATING FOR MULTINODE/REGION
        pth = fileparts(mfilename('fullpath'));
        pr = load([pth '/+fun/Priors2021a']);;
        
        else
            pr = load(priorfile);
        end
    
    % Extrinsics: restructure adjacency matrices
    %----------------------------------------------------------------------
    A     = DCM.A;
    A{1}  = A{1} | A{3};                              % forward
    A{2}  = A{2} | A{3};                              % backward

    % [F] SP -> SS & DP
    pE.A{1} = A{1}*32 - 32; % 0.5 = 1001.60944 - 1000
    pC.A{1} = A{1}/8;

    % [B] DP -> SP & SI
    pE.A{2} = A{2}*32 - 32; 
    pC.A{2} = A{2}/8;

    % [Extrinsic TP-TP] (none)
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
        
    % Average Firing
    pE.S = pr.pE.S;
    pC.S = pr.pC.S;
    
    % Average (Self) Excitation Per Population
    pE.G = zeros(ns,np);
    pC.G = zeros(ns,np); 
    
    % Average Membrane Capacitance
    pE.CV = pr.pE.CV;
    pC.CV = pr.pC.CV;
    
    % Average Background Activity
    %pE.E = 0;
    pE.E = pr.pE.E;
    pC.E = pr.pC.E;
    
    % Average Intrinsic Connectivity
    pE.H = repmat(pr.pE.H, [1 1 ns]);
    pC.H = repmat(pr.pC.H, [1 1 ns]);
             
    % Receptor Time Constants (Channel Open Times)
    pE.T = repmat(pr.pE.T,[ns,1]);
    pC.T = repmat(pr.pC.T,[ns,1]);
    
    % Parameters on input bump: delay, scale and width
    pE.R = pr.pE.R;
    pC.R = pr.pC.R;
    
    % Delays: states - intrinsic & extrinsic    
    pE.D = pr.pE.D;
    pC.D = pr.pC.D;
    
    %pE.D0 = [-0.0872 0.1470];
    pE.D0 = pr.pE.D0;
    pC.D0 = [1 1]*0;
    
    % Dipole position (not in use)
    pE.Lpos = sparse(3,0);
    pC.Lpos = sparse(3,0);
    
    % Electrode gain
    pE.L = repmat(pr.pE.L(1)         ,[ns,1]);
    pC.L = repmat(pr.pC.L(1)         ,[ns,1]);
        
    pE.J = pr.pE.J;
    pC.J = pr.pC.J;
        
    % Noise Components - a, b & c
    pE.a = pr.pE.a;
    pC.a = pr.pE.a;
    
    pE.b = pr.pE.b;
    pC.b = pr.pC.b;
    
    pE.c = repmat(pr.pE.c,[1 ns]);
    pC.c = repmat(pr.pC.c,[1 ns]);
          
    pE.d = pr.pE.d;
    pC.d = pr.pC.d;
       
    pE.h = repmat(pr.pE.h,  [ns 1]);
    pC.h = repmat(pr.pC.h,  [ns 1]);
    
    pE.m = repmat(pr.pE.m,  [ns 1]);
    pC.m = repmat(pr.pC.m,  [ns,1]);
    
    pE.TV = pr.pE.TV;
    pC.TV = pr.pC.TV;
    
    pE.Mh = pr.pE.Mh;
    pC.Mh = pr.pC.Mh;
    
    pE.Hh = pr.pE.Hh;
    pC.Hh = pr.pC.Hh;
    
    try
        pE.ID = pr.pE.ID;
        pC.ID = pr.pC.ID;
    catch
        pE.ID = zeros(1,8);
        pC.ID = zeros(1,8);
    end
    
    try
        pE.Hn = repmat(pr.pE.Hn, [1 1 ns]);
        pC.Hn = repmat(pr.pC.Hn, [1 1 ns]);
    catch
        pE.Hn = repmat(zeros(8,8), [1 1 ns]);
        pC.Hn = repmat(pC.H, [1 1 ns]);
    end
    
    try
        pE.pr = pr.pE.pr;
        pC.pr = pr.pC.pr;
    catch
        pE.pr = zeros(1,8);
        pC.pr = zeros(1,8);
    end
    
    % Pack pE & pC into M
    DCM.M.pE = pE;
    DCM.M.pC = pC;
    
end
    
    
end
