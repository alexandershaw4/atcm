function Qp = get_posteriors_tcm2024(Ep)

F = fixedvalues();

% intrinsics
Qp.AMPA = exp(Ep.H).*F.GEa;
Qp.GABAA = exp(Ep.H).*F.GIa;
Qp.NMDA = exp(Ep.Hn).*F.GEn;
Qp.GABAB = exp(Ep.Gb).*F.GIa;


% time constants in ms
Qp.KE  = 1./(exp(-Ep.T(:,1))*1000/2.2);
Qp.KI  = 1./(exp(-Ep.T(:,2))*1000/5);
Qp.KN  = 1./(exp(-Ep.T(:,3))*1000/100);
Qp.KB  = 1./(exp(-Ep.T(:,4))*1000/300);  
Qp.KM  = 1./(exp(-Ep.T(:,5))*1000/160);  
Qp.KH  = 1./(exp(-Ep.T(:,6))*1000/100);  

% thalamo-cortical and cortico-thalamic delays
Qp.CT = (8* exp(Ep.CT))/1000; %60;
Qp.TC = (3* exp(Ep.TC))/1000; %20;

% cell / axon conductance delays
ID = [2 1 1 1 1 2 1 2];
Qp.ID = (ID.*exp(Ep.ID)/1000); 

Qp.J = exp(Ep.J);

Qp.scale_NMDA = exp(Ep.scale_NMDA);

end


function F = fixedvalues()

F.SA   = [1   0   0   0   0;   %  SS    % added TP->SP
    0   1   0   0   0;   %  SP
    0   1   0   0   0;   %  SI
    0   0   0   0   0;   %  DP
    0   0   0   0   0;   %  DI
    0   0   0   0   0;   %  TP % 0 in ket study
    0   0   0   0   1;   %  rt % 0 in ket study
    0   0   0   1   0]/8;%  rc % 0 in ket study

F.SA(:,[3 4 5]) = 0; % For ket study

% % extrinsic NMDA-mediated connections (F B) - from superficial and deep pyramidal cells
% %--------------------------------------------------------------------------
F.SNMDA = [1   0   0   0   0;   %  SS
    0   1   0   0   0;   %  SP
    0   1   0   0   0;   %  SI
    0   0   0   0   0;   %  DP
    0   0   0   0   0;   %  DI
    0   0   0   0   0;   %  TP % 0 in ket study
    0   0   0   0   1;   %  rt % 0 in ket study
    0   0   0   1   0]/8;%  rc % 0 in ket study

F.SNMDA(:,[3 4 5]) = 0; % For ket study

% Excitatory (np x np): AMPA & NMDA
%--------------------------------------------------------------------------
F.GEa = [  0     0     0     0     0     2     0     2;
           2     0     0     0     0     0     0     0;
           0     2     0     0     0     0     0     0;
           0     2     0     0     0     0     0     0;
           0     0     0     2     0     0     0     0;
           0     0     0     2     0     0     0     0;
           0     0     0     0     0     0     0     2;
           2     0     0     0     0     2     0     0];

F.GEn =   [0     0     0     0     0     2     0     2;
           2     2     2     0     0     0     0     0;
           0     2     2     0     0     0     0     0;
           0     2     0     0     0     0     0     0;
           0     0     0     2     0     0     0     0;
           0     0     0     2     0     0     0     0;
           0     0     0     0     0     0     0     2;
           2     0     0     0     0     2     0     0];

% Inhibitory connections (np x np): GABA-A & GABA-B
%--------------------------------------------------------------------------
F.GIa =[8     0     10    0     0     0     0     0;
        0    18     10    0     0     0     0     0;
        0     0     10    0     0     0     0     0;
        0     0     0     8     6     0     0     0;
        0     0     0     0    14     0     0     0;
        0     0     0     0     6     8     0     0;
        0     0     0     0     0     0     8     0;
        0     0     0     0     0     0     8     8];

end


