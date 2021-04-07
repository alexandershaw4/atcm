function DCM = parameters6(DCM,Ns)

np = 6;

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

% % [Extrinsic T->T] rt -> rt & rl
% AA      = A{1} | A{2};
% 
% pE.A{4} = AA*32 - 32;
% pC.A{4} = AA/8;
% 
% % [Extrinsic T->T] rl -> rt & rl
% pE.A{5} = AA*32 - 32;
% pC.A{5} = AA/8;

% NMDA = AMPA
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

% exogenous inputs
%----------------------------------------------------------------------
C     = DCM.C;
C     = ~~C;
pE.C  = C*32 - 32;
pC.C  = C/8;

% new multi-input model
pE.C = repmat(pE.C,[1 3]);
pC.C = repmat(pC.C,[1 3]);


% Intrinsic (local) parameters, per region
%----------------------------------------------------------------------
ns = length(A{1}); % number of regions / nodes
np = 6;            % number of populations per region
nk = 7;            % number of states per population

% Average Firing
pE.S = zeros(1,np);
pC.S = ones(1,np)*0.0625;

% Average (Self) Excitation Per Population
pE.G = zeros(ns,np);
pC.G = zeros(ns,np);
pC.G(1:2) = 1;

% Average Membrane Capacitance
pE.CV =  [-0.0382646336383898
    0.00217382054723993
    0.0200769192441413
    0.0128455748594877
    0.0186176599910718
    -0.0369774537374279]';
pC.CV = ones(1,np) * 0.0625;

% Average Background Activity
pE.E = 0;
pC.E = 0;

% Average Intrinsic Connectivity
pE.H = repmat([...
    0.0073358603621668                         0        0.0100245782752032                         0                         0                         0                         0        -0.121912075196584;
    -0.023793524625747        0.0100227054481604        0.0283195175849927                         0                         0                         0                         0                         0;
    0        0.0207664577343551       -0.0204340813957294                         0                         0                         0                         0                         0;
    0                         0                         0       0.00933785756297514                         0                         0                         0                         0;
    0                         0                         0                         0                         0                         0                         0                         0;
    0                         0                         0       -0.0297913734339767                         0       -0.0219338394620789                         0      0.000587849049239577;
    0                         0                         0                         0                         0                         0                         0                         0;
    0                         0                         0                         0                         0      -0.00200828167128864                         0       -0.0109368241673351], [1 1 ns]);
%pC.H = ~~pE.H / 8;
pC.H =repmat([...
    0.001953125                         0               0.001953125                         0                         0                         0                         0               0.001953125;
    0.001953125               0.001953125               0.001953125               0.001953125                         0               0.001953125                         0               0.001953125;
    0               0.001953125                     0.002                         0                         0                         0                         0                         0;
    0                         0                         0               0.001953125                         0                         0                         0               0.001953125;
    0                         0                         0                         0                         0                         0                         0                         0;
    0                         0                         0               0.001953125                         0               0.001953125                         0               0.001953125;
    0                         0                         0                         0                         0                         0                         0                         0;
    0                         0                         0                         0                         0               0.001953125                         0               0.001953125], [1 1 ns]);

% CROP out thalamus
for i = 1:ns
    pEH(:,:,i) = (pE.H(1:6,1:6,i));
    pCH(:,:,i) = (pC.H(1:6,1:6,i));
end
pE.H = pEH;
pC.H = pCH;

pC.H(:,:,1) = eye(6)/8;
pC.H(:,:,2) = eye(6)/8;
pC.H(:,:,3) = eye(6)/8;


% Receptor Time Constants (Channel Open Times)
pE.T = repmat([0.00576645060668367 -0.0191870045945023 -0.00678195288168858 -0.0133612387162409],[ns,1]);
pC.T =  ( ones(ns,4) / 8);

% Parameters on input bump: delay, scale and width
pE.R = [0 0 0];
pC.R = [0 0 0];

% Delays: states - intrinsic & extrinsic
pE.D = [0 0];
pC.D = [0 0];

% Dipole position (not in use)
pE.Lpos = sparse(3,0);
pC.Lpos = sparse(3,0);

% Electrode gain
pE.L = repmat(4.97252413750747,[ns,1]);
pC.L = repmat(64              ,[ns,1]);

% Contributing states
J           = zeros(1,np,nk) - 1000;
J(:,[1 2 4 6],1)    = log([.2 .8 .2 .2]);
J(isinf(J)) = -1000;
pE.J        = spm_vec(J)';
pC.J        = spm_vec(J)' * 0;

% Noise Components - a, b & c
pE.a = repmat([-0.0265936894424734;0],[1 ns]);
pC.a = repmat([0.03125;0]            ,[1 ns]);

pE.b = [0;0];
pC.b = [0;0];

pE.c = repmat([-4.99912119824703;-0.0902806280629704],[1 ns]);
pC.c = repmat([0.03125;0.03125],                      [1 ns]);

% Neuronal innovations
pE.d = [ -0.748515803154631;
    -0.32731663446822;
    0.0101072993361102;
    -0.242119529214598;
    -0.00320081426472445;
    -0.00550311277258295;
    -0.000392685304729754;
    0.00472838376680181]';

pC.d = ones(length(pE.d),1) / 8;

pE.h = repmat(0,  [ns 1]);
pC.h = repmat(0,  [ns 1]);

pE.m = repmat(0,  [ns 1]);
pC.m = repmat(0,  [ns,1]);

pE.TV = [0.0166517239998525        -0.046353787602425        0.0107371114589147      -0.00533791556468411     -0.000487422444009102         0.133214597172099];
pC.TV = ones(1,np)*0.125;

pE.Mh = zeros(np,1);
pC.Mh = zeros(np,1);

pE.Hh = [0 0];
pC.Hh = [0 0];

% Pack pE & pC into M
DCM.M.pE = pE;
DCM.M.pC = pC;
