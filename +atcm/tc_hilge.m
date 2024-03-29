function [f,J] = tc_hilge(x,u,P,M,fso)
% State equations for an extended canonical thalamo-cortical neural-mass model.
%
% This model implements a conductance-based canonical thalamo-cortical circuit,
% with cytoarchitecture inspired by Gilbert & Wiesel (1983), Douglas & 
% Martin (2004) and Traub (2004) models.
%
% The equations of motion are Moris Lecar-esque equations, similar to Moran
% (2011), but with conductances for AMPA, NMDA, GABA-A, & GABA-B channels. 
% These 'channels' feature their own reversal poentials and rate constants:
%
% K  = -70           (Leak)
% Na =  60  / 4 ms   (AMPA)
% Cl = -90  / 16 ms  (GABA-A)
% Ca =  10  / 100 ms (NMDA)   + voltage mag switch
% B  = -100 / 200 ms (GABA-B)
% f  = -40
%
% FORMAT [f,J,Q,D] = atcm.tcm_hilge(x,u,P,M)
%
% x - states and covariances
%
% x(i,j,k)        - k-th state of j-th population of i-th source
%                   i.e., running over sources, pop. and states
%
%   population: 1  - Spint stellates (L4)
%               2  - Superficial pyramids (L2/3)
%               3  - Inhibitory interneurons (L2/3)     
%               4  - Deep pyramidal cells (L5)
%               5  - Deep interneurons (L5)
%               6  - Thalamic projection neurons (pyramid) (L6)
%               7  - Reticular cells (Thal)
%               8  - Thalamo-cortical relay cells (Thal)
%
%
%        state: 1 V   - voltage
%               2 gE  - conductance: AMPA   (excitatory)
%               3 gI  - conductance: GABA-A (inhibitory)
%               4 gN  - conductance: NMDA   (excitatory)
%               5 gB  - conductance: GABA-B (inhibitory)
%               6 gM  - conductance: M-channels (inhibitory)
%               7 gih - conductance: H-channels (inhibitory)
%
%      outputs: f = model states as a vector - hint: spm_unvec(f,M.x) 
%               J = system Jacobian - dfdx
%               Q = delay operator  - Q = inv(1 - D.*dfdx)*f(x(t))
%               D = states delay matrix
%
% Info:
%  - Ih is a hyperpolarization-activated cation current mediated by HCN channel
%  - non-selective, voltag gated, responsible for cariac 'funny' (pacemaker) current
%
%  - M-channels (aka Kv7) are noninactivating potassium channels
%  - M is unique because it is open at rest and even more likely to be open during depolarization
%  - M is a pip2 regulated ion channel
%
% Alexander Shaw 2019: ShawA10@cardiff.ac.uk
%
% Notes, changes, updates:
%
% Extrinsics connection matrices [ampa but AN{n} is nmda equiv]:
% A{1} = Forward  SP -> SS & DP
% A{2} = Backward DP -> SP & SI
% A{3} = Back/Lat TP -> SS & TP
% A{4} = Inter-Thal [B] RT -> RC
% A{5} = Inter-Thal [F] RC -> RT
%
% Dr Alexander Shaw | 2020 | alexandershaw4[@]gmail.com


% Flag: include M- & H- channels on L6 TP & Thalamic Relay cells, or not
%--------------------------------------------------------------------------
IncludeMH = 1;


inputu = u;

% if isfield(P,'inputs') && length(u) == 1
%     u = repmat(u,[8 1]);
% else u = u(:);
% end
 
% get dimensions and configure state variables
%--------------------------------------------------------------------------
ns   = size(M.x,1);                      % number of sources
np   = size(M.x,2);                      % number of populations per source
nk   = size(M.x,3);                      % number of states per population
x    = reshape(x,ns,np,nk);              % hidden states 


% extrinsic connection strengths
%==========================================================================
 
% exponential transform to ensure positivity constraints
%--------------------------------------------------------------------------
for i = 1:length( P.A )
    A{i}  = exp(P.A{i});
    AN{i} = exp(P.AN{i});
end

C     = exp(P.C); 
 

% detect and reduce the strength of reciprocal (lateral) connections
%--------------------------------------------------------------------------
for i = 1:length(A)
    L    = (A{i} > exp(-8)) & (A{i}' > exp(-8));
    A{i} = A{i}./(1 + 8*L);
end

            
% intrinsic connection strengths
%==========================================================================
G    = full(P.H);
G    = exp(G);

% nmda matrix
%P.Hn=P.H;
%P.Hn(3,3)=P.nmda;

Gn = full(P.Hn);
Gn = exp(Gn);
%Gn = exp(full(G));

% this was for specifying trial-specific intrinsic connections (w Betas)
% in the LTP project (Sumner, Spriggs, Shaw 2020)
%--------------------------------------------------
% if ~all(size(P.G)==np)
%     for i = 1:size(G,3)
%         % trial specific intrinsic effects !
%         Gtrial   = diag( (P.G));
%         G(:,:,i) = G(:,:,i) + Gtrial; 
%     end
% elseif all(size(P.G)==np)
%     % a full 8*8 connectivity for this trial
%     for i = 1:size(G,3)
%         G(:,:,i) = G(:,:,i) + P.G;
%     end
% end




% connectivity switches
%==========================================================================
% 1 - excitatory spiny stellate cells (granular input cells)
% 2 - superficial pyramidal cells     (forward  output cells)
% 3 - inhibitory interneurons         (intrisic interneuons)
% 4 - deep pyramidal cells            (backward output cells)
% 5 - deep interneurons               
% 6 - thalamic projection pyramidal cells (with m- and h- currents)
% 7 - thalamic reticular cells
% 8 - thalamic relay cells (with m- and h- currents)
%
% Thalamic cells attached to different cortical regions (models) are laterally connected


% % extrinsic connections (F B) - from superficial and deep pyramidal cells
% %--------------------------------------------------------------------------
%       SP  DP  tp  rt  rc
SA   = [1   0   0   0   0;   %  SS    % added TP->SP
        0   1   0   0   0;   %  SP
        0   1   0   0   0;   %  SI
        0   0   0   0   0;   %  DP
        0   0   0   0   0;   %  DI
        0   0   0   0   0;   %  TP % 0 in ket study
        0   0   0   0   1;   %  rt % 0 in ket study
        0   0   0   1   0]/8;%  rc % 0 in ket study
    
    SA(:,[3 4 5]) = 0; % For ket study
    
    
% % extrinsic NMDA-mediated connections (F B) - from superficial and deep pyramidal cells
% %--------------------------------------------------------------------------    
SNMDA = [1   0   0   0   0;   %  SS
         0   1   0   0   0;   %  SP
         0   1   0   0   0;   %  SI
         0   0   0   0   0;   %  DP
         0   0   0   0   0;   %  DI
         0   0   0   0   0;   %  TP % 0 in ket study
         0   0   0   0   1;   %  rt % 0 in ket study
         0   0   0   1   0]/8;%  rc % 0 in ket study

     SNMDA(:,[3 4 5]) = 0; % For ket study
     
% intrinsic connectivity switches
%--------------------------------------------------------------------------    
%   population: 1  - Spint stellates (L4)                : e
%               2  - Superficial pyramids (L2/3)         : e
%               3  - Inhibitory interneurons (L2/3)      : i
%               4  - Deep pyramidal cells (L5)           : e
%               5  - Deep interneurons (L5)              : i
%               6  - Thalamic projection neurons -L6     : e
%               7  - Reticular cells (Thal)              : i
%               8  - Thalamo-cortical relay cells (Thal) : e

GEa = zeros(8,8);
GIa = zeros(8,8);

% Excitatory (np x np): AMPA & NMDA
%--------------------------------------------------------------------------


GEa = [  0     0     0     0     0     2     0     2;
         4     4     0     0     0     0     0     0;
         0     2     0     0     0     0     0     0;
         0     2     0     0     0     0     0     0;
         0     0     0     2     0     0     0     0;
         0     0     0     8     0     0     0     0;
         0     0     0     0     0     0     0     2;
         2     0     0     0     0     2     0     0];

% GEa = [  0     0     0     0     0     4     0     2;
%          4     4     0     0     0     4     0     0;
%          0     2     0     0     0     0     0     0;
%          0     2     0     0     0     2     0     0;
%          0     0     0     2     0     0     0     0;
%          2     2     0     0     0     0     0     0;
%          0     0     0     0     0     0     0     2;
%          2     0     0     4     0     2     0     0];

GEa = GEa*2;

%GEa = GEa/8;


GEn = GEa;
%GEn = GEn + eye(8)/8;



% Inhibitory connections (np x np): GABA-A & GABA-B
%--------------------------------------------------------------------------

si = 8;
GIa =  [ si    0     8     0     0     0     0     0;
         0     si    64    0     0     0     0     0;
         0     0     si    0     0     0     0     0;
         0     0     0     12    4     0     0     0;
         0     0     16    0     4     0     0     0;
         0     0     0     0     8     4     0     0;
         0     0     0     0     0     0     12    0;
         0     0     0     0     0     0     4     si];
     
GIa=GIa*2;

%GIa = GIa/8;

% GIbg = ~GIa/32;
% GEbg = ~GEa/32;
% 
% %GIbg(:,[1 2 4 6 8])=0;
% %GEbg(:,[3 5 7]) = 0;
% 
% GIa = GIa/2 + ( GIbg * exp(P.TV(1)));
% GEa = GEa/2 + ( GEbg * exp(P.TV(2)));


GIb = GIa;

if isfield(P,'scale')
    GEa = GEa * exp(P.scale(1));
    GEn = GEn * exp(P.scale(2));
    GIa = GIa * exp(P.scale(3));
    GIb = GIb * exp(P.scale(4));
end


if IncludeMH
    
    % M- & H- channel conductances (np x np) {L6 & Thal Relay cells only}
    %----------------------------------------------------------------------
    % https://www.sciencedirect.com/science/article/pii/S0006349599769250
    VM   = -70;                            % reversal potential m-channels          
    VH   = -30;                            % reversal potential h-channels 

%     GIm  = sparse([6 8],[6 8],4,8,8);
%     %GIm  = eye(8)*4;
%     Mh   = diag(exp(P.Mh));
% 
%     GIh      = full(sparse([6 8],[6 8],4   ,8,8)); % 1/4
%     Hh       = exp(P.Hh);
%     GIh(6,6) = GIh(6,6)*Hh(1);
%     GIh(8,8) = GIh(8,8)*Hh(2);


    GIm = diag(4*[0 0 0 0 0 1 0 1].*exp(P.Mh(:)'));
    GIh = diag(4*[0 0 0 0 0 1 0 1].*exp(P.Hh(:)'));

    KM    = (exp(-P.T(:,5))*1000/160) ;               % m-current opening + CV
    KH    = (exp(-P.T(:,6))*1000/100) ;               % h-current opening + CV
    h     = 1 - spm_Ncdf_jdw(x(:,:,1),-100,300); % mean firing for h-currents
    %h      = 1./(1+exp((x(:,:,1)+81)/7));
end

% Channel rate constants [decay times]
%--------------------------------------------------------------------------
KE  = exp(-P.T(:,1))*1000/2.2;%3;            % excitatory rate constants (AMPA) % 2 to 5
KI  = exp(-P.T(:,2))*1000/5;%6;           % inhibitory rate constants (GABAa)
KN  = exp(-P.T(:,3))*1000/100;%40;          % excitatory rate constants (NMDA) 40-100
KB  = exp(-P.T(:,4))*1000/300;          % excitatory rate constants (NMDA)

% cojuld even use number from this friston paper
%https://www.sciencedirect.com/science/article/pii/S0361923000004366?via%3Dihub
% ampa = 1.2 to 2.4 ms
% gabaa -   6ms
% nmda - 50 ms
%KN  = exp(-P.T(:,3))*1000/50;    

% gaba-b maybe evern 300 or 500ms

% now using faster AMPA and GABA-A dynamics based on this book:
% https://neuronaldynamics.epfl.ch/online/Ch3.S1.html#:~:text=GABAA%20synapses%20have%20a,been%20deemed%203%20times%20larger.

%KE  = exp(-P.T(:,1))*1000/3;            % excitatory rate constants (AMPA)
%KN  = exp(-P.T(:,3))*1000/150;          % excitatory rate constants (NMDA)
%KI  = exp(-P.T(:,2))*1000/6;           % inhibitory rate constants (GABAa)

%KE = KE * exp( P.TT(1) );
%KI = KI * exp( P.TT(2) );
%KN = KN * exp( P.TT(3) );
%KB = KB * exp( P.TT(4) );


% Trial effects on time constants: AMPA & NMDA only
if isfield(P,'T1')
    KE = KE + P.T1(1);
    KN = KN + P.T1(2);
end

% Voltages [reversal potentials] (mV)
%--------------------------------------------------------------------------
VL   = -70;                               % reversal  potential leak (K)
VE   =  60;                               % reversal  potential excite (Na)
VI   = -90;                               % reversal  potential inhib (Cl)
VR   = -52;   %55                            % threshold potential (firing)
VN   =  10;                               % reversal Ca(NMDA)   
VB   = -100;                              % reversal of GABA-B

%VE = VE * exp(P.pr(1));
%VI = VI * exp(P.pr(2));
%VN = VN * exp(P.pr(3));
%VB = VB * exp(P.pr(4));

VR = VR * exp(P.pr(1));

% membrane capacitances {ss  sp  ii  dp  di  tp   rt  rl}
%--------------------------------------------------------------------------
CV   = exp(P.CV).*      [128*3 128 128/2 128 64  128  64  64*2]/1000;  

% leak conductance - fixed
%--------------------------------------------------------------------------
GL   = 1 ;       

% mean-field effects:
%==========================================================================

% neural-mass approximation to covariance of states: trial specific
%----------------------------------------------------------------------
%Vx   = exp(P.S)*32; % 32
%    {ss    sp    ii    dp    di    tp    rt    rl}
%Vx = [Vx(1) Vx(1) Vx(2) Vx(1) Vx(2) Vx(1) Vx(3) Vx(3)];

    %-Approximate integral
    %--------------------------------------------------------------------------
    %x0    = (x(:,:,1) - VR)./sqrt(abs(Vx));
    %F    = sqrt(1 - exp(-(2/pi)*x0.^2))/2;
    %i    = x0 < 0;
    %F(i) = -F(i);
    %m    = F + 1/2;
    
    %SPf = 1./(1 + exp(-exp(P.SP).*(x(:,2,1)-VR)));
    %m(2)  = SPf';
        
    FF = 1./(1 + exp(-exp(P.S).*(x(:,:,1)-VR)));
    
    % probability of firing: 
    % [-70 : -55] mV == LHS of Gaussian co VR
    % [-55 : +30] mV == 1
    % [ >= 30 ]   mV == 0
        
    % do this to map the probability of firing over V:
    % figure,plot(-70:1:70,pdf(makedist('normal',VR,2),(-70:1:70)) ./ pdf(makedist('normal',VR,2),VR))
    
     % reset to refractory
    %RS = find( x(:,:,1) >= 30 );
    %if any(RS)
    %    x(1,RS,1) = -90;
    %end
        
    %FF =  pdf(makedist('normal',VR,2),x(:,:,1)) ./ pdf(makedist('normal',VR,2),VR);
    
    %VVR = VR*exp(P.S); % try 'Lognormal'
    %for i = 1:8
    %    FF(i) =  pdf(makedist('normal',VVR(i),2),x(:,i,1)) ./ pdf(makedist('normal',VVR(i),2),VVR(i));
    %end
    
    %x(1,:,1) = squeeze(x(1,:,1)) .* (1./(1 + exp( - exp(P.S(:)) )))';
        
    RS = 30 ;
    
    Fu = find( x(:,:,1) >= VR );
    FF(Fu) = 1;
    
    Fl = find( x(:,:,1) >= RS );
    FF(Fl) = 0;

    
    m  = FF;%.*exp(P.S);

    
% extrinsic effects
%--------------------------------------------------------------------------
a       = zeros(ns,5);
an      = zeros(ns,5); 
a(:,1)  = A{1}*m(:,2);                      % forward afference  AMPA - SP->SS
a(:,2)  = A{2}*m(:,4);                      % backward afference AMPA - DP->SP
a(:,3)  = A{3}*m(:,6);                      % FWD thalamic projection pyramids
a(:,4)  = A{4}*m(:,7);                      % LAT reticular AMPA
a(:,5)  = A{5}*m(:,8);                      % LAT relay AMPA
an(:,1) = AN{1}*m(:,2);                     % forward afference  NMDA
an(:,2) = AN{2}*m(:,4);                     % backward afference NMDA
an(:,3) = AN{3}*m(:,6);                     % thalamic projection pyramids
an(:,4) = AN{4}*m(:,7);                     % reticular NMDA
an(:,5) = AN{5}*m(:,8);                     % relay NMDA

% Averge background activity and exogenous input
%==========================================================================
BE     = exp(P.E)*0.8;

% input(s)
%--------------------------------------------------------------------------
% if isfield(M,'u')
%       U =   u;%(:); % endogenous input
% else; U = C*u;%(:); % exogenous input
% end

% flow over every (ns x np) subpopulation
%==========================================================================
f     = x;

% Thalamo-cortical flow [eq. motion] over modes, populations, states...
%--------------------------------------------------------------------------
for i = 1:ns
   
        % allow switching of thalamus on/off
        %HasThal = M.HasThal(i);
            
        % input scaling: 
        %------------------------------------------------------------------
        %if any(full(U(:))) ;
            dU = u(:)*C(i,1);
        %else
        %    dU = 0;
        %end
                
        Gsc = ~eye(8);
        Gsc = Gsc +  (diag(exp(P.Gsc)));
                
        % CT and TC delays on G
        %-----------------------------------------------------------------
         CT = 8*exp(P.CT); %60;
         TC = 3*exp(P.TC); %20;

        %CT = 0.06*exp(P.CT); %60;
        %TC = 0.02*exp(P.TC); %20;

        tcd = zeros(8,8);

        tcd(1:6,[7 8]) = TC;
        tcd([7 8],1:6) = CT;

        % apply operator to G and GN
        G  = (1./(1+tcd)).*G;
        Gn = (1./(1+tcd)).*Gn;
                
        % intrinsic coupling - parameterised
        %------------------------------------------------------------------
        E      = ( (G(:,:,i).*GEa).*Gsc )*m(i,:)'; % AMPA currents
        ENMDA  = (Gn(:,:,i).*GEn)*m(i,:)'; % NMDA currents
        I      = ( G(:,:,i).*GIa)*m(i,:)'; % GABA-A currents
        IB     = ( G(:,:,i).*GIb)*m(i,:)'; % GABA-B currents
                
        if IncludeMH
            
            % intrinsic coupling - non-parameterised: intrinsic dynamics
            %--------------------------------------------------------------
            Im     = GIm*m(i,:)'; % M currents
            Ih     = GIh*h(i,:)'; % H currents
        end
        
        % extrinsic coupling (excitatory only) and background activity
        %------------------------------------------------------------------
        E     = (E     +  BE  + SA   *a (i,:)')*2;
        ENMDA = (ENMDA +  BE  + SNMDA*an(i,:)')*2;
                
        
        % and exogenous input(U): 
        %------------------------------------------------------------------
        %if length(u) == 1
            
        %input_cell        = [8];
        %E(input_cell)     = E(input_cell) + dU;

        % flag for the oscillation injection desribed here: 
        % https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5310631/
        
        if length(u) > 1
            E(8) = E(8) + dU(1);
            E(6) = E(6) + dU(2);
        else
            input_cell        = [8 1];
            E(input_cell)     = E(input_cell) + dU;
        end
       
        %E(2) = E(2) + exp(P.a(1)) * m(2);
        %I(3) = I(3) + exp(P.a(2)) * m(3);
        %E(4) = E(4) + exp(P.a(3)) * m(4);
        %E(1) = E(1) + exp(P.a(4)) * m(1);


         %   ENMDA(input_cell)     = ENMDA(input_cell) + dU*10;
        %elseif length(u) == 56
            %E(1:8) = E(1:8) + u(1:8)*10*C(1);
            %ENMDA(1:8) = ENMDA(1:8) + u(1:8)*10*C(1);

        %    E([1 2 4 6 8]) = E([1 2 4 6 8]) + (u([1 2 4 6 8])*10*C(1));
        %    ENMDA([1 2 4 6 8]) = ENMDA([1 2 4 6 8]) + (u([1 2 4 6 8])*10*C(2));
        %    I([3 5 7]) = I([3 5 7]) + (u([3 5 7])*10*C(3)); 
            

        %end
       % ENMDA(input_cell) = ENMDA(input_cell) * dU';
        
        %I(7) = I(7) * dU;
        %E(7) = E(7) - I(7);

%         if nargin > 4 && ~isempty(fso)
%             I(3) = I(3) + fso;
%             E(3) = E(3) + fso;
% 
%            % I(2) = I(2) + fso;
%         end
                       
        % Voltage equation
        %==================================================================
        if ~IncludeMH
            
          f(i,:,1) =         (GL*(VL - x(i,:,1))+...
                       1.0*x(i,:,2).*(VE - x(i,:,1))+...
                       1.0*x(i,:,3).*(VI - x(i,:,1))+...
                       1.0*x(i,:,5).*(VB - x(i,:,1))+...
                       1.0*x(i,:,4).*(VN - x(i,:,1)).*mg_switch(x(i,:,1)))./CV;
            
        elseif IncludeMH
            
          
          f(i,:,1) =  (GL*(VL - x(i,:,1))+...
                       x(i,:,2).*((VE - x(i,:,1)))+...
                       x(i,:,3).*((VI - x(i,:,1)))+...
                       x(i,:,5).*((VB - x(i,:,1)))+...
                       x(i,:,6).*((VM - x(i,:,1)))+...
                       x(i,:,7).*((VH - x(i,:,1)))+...
                       x(i,:,4).*((VN - x(i,:,1))).*mg_switch(x(i,:,1)))./CV;
        end
                   
        % Conductance equations
        %==================================================================           
        
        f(i,:,2) = (E'     - x(i,:,2)).* (KE(:,:)');%*pop_rates);
        f(i,:,3) = (I'     - x(i,:,3)).* (KI(:,:)');%*gabaa_rate);
        f(i,:,5) = (IB'    - x(i,:,5)).* (KB(:,:)');%*pop_rates);
        f(i,:,4) = (ENMDA' - x(i,:,4)).* (KN(:,:)');%*nmdat);
        
        if IncludeMH
            f(i,:,6) = (Im'    - x(i,:,6)).*(KM(i,:) );%*pop_rates );
            f(i,:,7) = (Ih'    - x(i,:,7)).*(KH(i,:) );%*pop_rates );
        end
                
        %f(i,:,7) = FF;
        
        % Conductance Delays: f = state update, x = state previous
        % population delays
        %------------------------------------------------------------------
        % short range neurons are around 0.014m length and conduct around
        % 120 m/s so ~1.68
        
%         df = f(:) - x(:);
%         d  = exp(P.ID).*[1 1 1 1 1 1 1 1];
%         %d=[1 1 1 1 1 1 1 1];
%         s = 1;
%         l = 1.68;
%         
%         d = [l l s l s s s s]./d;
%         % exp(-P.T(:,1))*1000/2.2;
% 
%         %d = repmat(d(:),[nk 1]);
%         
% 
%         f = spm_unvec( spm_vec(x) + df(:).*repmat(d(:),[nk 1]), f);

                                
%         % receptor and cell specific delays
%         df = f - x;
%         
%         % SP AMPA
%         f(:,2,2) = x(:,2,2) + df(:,2,2) * exp(P.delays(1));
%         
%         % SP NMDA
%         f(:,2,4) = x(:,2,4) + df(:,2,4) * exp(P.delays(2));
%         
%         % SP GABA
%         f(:,2,3) = x(:,2,3) + df(:,2,3) * exp(P.delays(3));
%         
%         % SI GABA
%         f(:,3,3) = x(:,3,3) + df(:,3,3) * exp(P.delays(4));
%         
%         % DP AMPA
%         f(:,4,2) = x(:,4,2) + df(:,4,2) * exp(P.delays(5));
%         
%         % SI NMDA
%         f(:,3,4) = x(:,3,4) + df(:,3,4) * exp(P.delays(6));
        

    
    

end


% vectorise equations of motion
%==========================================================================
f = spm_vec((f));
pE = P;
 
[J,Q,D]=deal([]);

if nargout < 2 || nargout == 5, return, end

% Only compute Jacobian (gradients) if requested
%==========================================================================
J = spm_cat(spm_diff(M.f,x,u,P,M,1));



if nargout < 3, return, end

% Only compute Delays if requested
%==========================================================================
% Delay differential equations can be integrated efficiently (but 
% approximately) by absorbing the delay operator into the Jacobian
%
%    dx(t)/dt     = f(x(t - d))
%                 = Q(d)f(x(t))
%
%    J(d)         = Q(d)df/dx
%--------------------------------------------------------------------------
% [specified] fixed parameters
%--------------------------------------------------------------------------
D  = [1 16];
d  = -D.*full(exp(P.D(1:2)))/1000;
Sp = kron(ones(nk,nk),kron( eye(np,np),eye(ns,ns)));  % states: same pop.
Ss = kron(ones(nk,nk),kron(ones(np,np),eye(ns,ns)));  % states: same source

% Thalamo cortical interactions: ~80ms round trip: 20 ms T->C, 60 ms C->T
%--------------------------------------------------------------------------
%Thalamocortical connections and forward connections from Vp to Vs had
%a mean delay of 3 ms, while corticothalamic connections and backward
%connections from Vs to Vp had a mean delay of 8 m - Neural Dynamics in a Model of the
%Thalamocortical System. I. Layers, Loops and the Emergence of Fast Synchronous Rhythms
% Lumer et al 1997

CT = 8; %60;
TC = 3; %20;

%CT = 60;
%TC = 20;

Tc              = zeros(np,np);
Tc([7 8],[1:6]) = CT  * exp(P.D0(1)); % L6->thal
Tc([1:6],[7 8]) = TC  * exp(P.D0(2)); % thal->ss
%Tc = Tc.*~~(GEa | GIa);
Tc = -Tc / 1000;
Tc = kron(ones(nk,nk),kron(Tc,eye(ns,ns)));


% ID = [4 1/4 1 8 1/2 4 2 20]/8;%2.4;
% ID = [1 1 1 1 1 1 1 20];
% ID = -ID.*exp(P.ID)/1000;
% 
% ID = repmat(ID,[1 nk]);
% 
% 
% Tc = Tc + diag(ID);


% Mean intra-population delays, inc. axonal etc. Seem to help oscillation
%--------------------------------------------------------------------------
Dp = ~Ss;                            % states: different sources
Ds = ~Sp & Ss;                       % states: same source different pop.
%Ds = Ds.*(~(Ds & Tc));              % remove t-c and c-t from intrinsic

if ~isfield(P,'delays')
    D  = d(2)*Dp + d(1)*Ds ;%+ Tc  ;
else
    D = d(1)*Ds + Tc  ;       %+ Dself;% Complete delay matrix
end

%D = d(2)*Dp + Tc; %%%%%!!!!!!

% Implement: dx(t)/dt = f(x(t - d)) = inv(1 - D.*dfdx)*f(x(t))
%                     = Q*f = Q*J*x(t)
%--------------------------------------------------------------------------
Q  = spm_inv(speye(length(J)) - D.*J);
%Q  = spm_inv(D.*J);
