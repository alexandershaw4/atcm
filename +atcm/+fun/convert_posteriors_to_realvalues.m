function p=convert_posteriors_to_realvalues(p)


% AMPA / NMDA connection switches
GEa(1,:) = [0   0   0   0   0   2   0   2]/1;
GEa(2,:) = [4   0   0   0   0   0   0   1]/1;
GEa(3,:) = [4   4   0   0   0   0   0   0]/1; 
GEa(4,:) = [0   4   0   0   0   2   0   1]/1;
GEa(5,:) = [0   0   0   4   0   2   0   0]/1;
GEa(6,:) = [0   0   0   2   0   0   0   1/4]/1; 
GEa(7,:) = [0   0   0   1   0   2   0   2]/1; 
GEa(8,:) = [0   0   0   1   0   2   0   0]/1;
GEa = GEa/8;
GEa = GEa .* ~eye(8);
GEa = GEa + eye(8);      % KILLED
GEn = GEa;
GEn = GEn + (eye(8)/8);

% GABA A/B connection switches
GIa(1,:) = [8   0   8   0   0   0   0   0 ];
GIa(2,:) = [0   16  32  0   0   0   0   0 ];
GIa(3,:) = [0   0   32  0   32  0   0   0 ];
GIa(4,:) = [0   0   0   8   12  0   0   0 ];
GIa(5,:) = [0   0   32  0   16  0   0   0 ];
GIa(6,:) = [0   0   0   0   32  8   0   0 ];
GIa(7,:) = [0   0   0   0   0   0   32  0 ];
GIa(8,:) = [0   0   0   0   0   0   8   32]; % 32!! SEPT
GIa = GIa/2;
GIb = GIa;

% convert - exp(p)*x
p.H = exp(p.H).*(GEa+GIa);

% NMDA
p.Hn = exp(p.Hn).*GEn;

% time constants - convert to ms
p.T(1)  = 1./(exp(-p.T(1))*1000/4);            % excitatory rate constants (AMPA)
p.T(2)  = 1./(exp(-p.T(2))*1000/16);           % inhibitory rate constants (GABAa)
p.T(3)  = 1./(exp(-p.T(3))*1000/100);          % excitatory rate constants (NMDA)
p.T(4)  = 1./(exp(-p.T(4))*1000/200);          % excitatory rate constants (NMDA)
p.m  = (exp(-p.m)*1000/160) ;               % m-current opening + CV
p.h  = (exp(-p.h)*1000/100) ;               % h-current opening + CV

% memebrane capacitance
p.CV = exp(p.CV).*[128 128 64  128 64  128  64  64*2]/1000;  

% population variances on synaptic delays
pop_rates = [1 1 1 1 1 1 1 1];
p.pr = pop_rates.*exp(p.pr);

% thalamic delays
CT = 8; 
TC = 3; 
p.D0(1) = CT  * exp(p.D0(1)); % L6->thal
p.D0(2) = TC  * exp(p.D0(2)); % thal->ss

% intrinsic delays
p.ID = exp(p.ID).*[2 1/4 1/2 4 1/2 4 2 2]/2;


        
        
        
        
        
        
        
        
        