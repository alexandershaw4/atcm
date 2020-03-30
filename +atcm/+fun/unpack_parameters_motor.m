function p  = unpack_parameters_motor(DCM,pE)
% returns parameters of the thalamo-cortical model in real-world units -
% i.e. seconds, mA, siemens etc.
%
%
%
% AS


% Model state space dimensions
[ns,npp,nk] = size(DCM.M.x);

% extrinsic adjacency matrices, trial specific connections and input
% strengths just need exponentiating
%--------------------------------------------------------------------------
p.A  = spm_unvec( exp(spm_vec(pE.A )) , pE.A  );
p.AN = spm_unvec( exp(spm_vec(pE.AN)) , pE.AN );
p.B  = spm_unvec(    (spm_vec(pE.B )) , pE.B  );
p.BN = spm_unvec(    (spm_vec(pE.BN)) , pE.BN );
p.C  = spm_unvec( exp(spm_vec(pE.C ))'.*[1 1/16 1/16 1/16] , pE.C  );

% intrisic connection parameters need exp & multiplying by their fixed vals
%--------------------------------------------------------------------------
H     = exp(pE.H);
G     = exp(pE.G);
E(1,:) = [0   1   0   0   0   4   0   2]/1;
E(2,:) = [4   0   0   1   0   0   0   0]/1;
E(3,:) = [4   4   0   0   0   0   0   0]/1; 
E(4,:) = [0   4   0   0   0   0   0   2]/1;
E(5,:) = [0   0   0   4   0   0   0   0]/1;
E(6,:) = [0   0   0   4   0   0   0   0]/1;
E(7,:) = [0   0   0   0   0   0   0   2]/1; 
E(8,:) = [0   0   0   0   0   2   0   0]/1;

I(1,:) = [8   0   2   0   0   0   0   0 ];
I(2,:) = [0   16  16  0   0   0   0   0 ];
I(3,:) = [0   0   32  0   0   0   0   0 ];
I(4,:) = [0   0   0   8   8   0   0   0 ];
I(5,:) = [0   0   0   0   16  0   0   0 ];
I(6,:) = [0   0   0   0   8   8   0   0 ];
I(7,:) = [0   0   0   0   0   0   32  0 ];
I(8,:) = [0   0   0   0   0   0   32   32];

for i  = 1:ns
    
    % intrinsc connections and inhibitory self-mods
    p.Intrinsic(:,:,i) = ( E.*H(:,:,i) ) + ( I.*H(:,:,i) );
    
end

% receptor opening and population conduction velocity time constants
%--------------------------------------------------------------------------
KE  = exp(-pE.T(:,1))*1000/4;            % excitatory rate constants (AMPA)
KI  = exp(-pE.T(:,2))*1000/16;           % inhibitory rate constants (GABAa)
KN  = exp(-pE.T(:,3))*1000/100;          % excitatory rate constants (NMDA)
KB  = exp(-pE.T(:,4))*1000/200;          % excitatory rate constants (NMDA)

p.K_AMPA = 1./KE';
p.K_GABA = 1./KI';
p.K_NMDA = 1./KN';
p.K_GABB = 1./KB';

%p.K_Mcur =  1./( exp(-pE.m)*1000/160 );   % m-current rate constant
%p.K_Hcur =  1./( exp(-pE.h)*1000/100 );   % h-current rate constant

% Temporal slowing
%DV         = 1./[2 1 1 2.2 1 2 1 2];
DV         = 1./[1 1 1 1   1 1 1 1];
p.DV       = DV.*exp(pE.TV);
        
% Membrane capacitance, Firing [S] & Background [BE]
%--------------------------------------------------------------------------
%p.CV   = exp(pE.CV).* [128 32  32  128 64  128  256 64*8]/1000;  
p.CV   = exp(pE.CV).* [128 32  32  128 64  128  256 64*2]/1000;  

p.S    = exp(pE.S)*32;
p.BE   = exp(pE.E)*0.8;

% Copy over custom/optional trial specific matrices
%--------------------------------------------------------------------------
if isfield(pE,'T1')
    p.T1 = pE.T1;
end

if isfield(pE,'G')
    p.G = pE.G;
end

if isfield(pE,'D0')
%     d0(1) = 60  * exp(pE.D0(1)); % L6->thal
%     d0(2) = 20  * exp(pE.D0(2)); % thal->ss
%     p.D0  = -d0 / 1000;
    
    %p.D0 = [60 20].*exp( pE.D0 );
    p.D0 = [8 3].*exp( pE.D0 );

end


    