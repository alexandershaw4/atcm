function components = iunpack(DCM,pE)
% integrate and unpack the returned components into a data structure
%
% components = iunpack(DCM,pE)
%

if nargin < 2 || isempty(pE)
    pE = DCM.M.pE;
end

% integrate
[y,w,s,g,t,pst,l,n,f] = feval(DCM.M.IS,pE,DCM.M,DCM.xU);

c.cells = {'ss' 'sp' 'si' 'dp' 'di' 'tp' 'rt' 'rl'};

% precompute function connectivity for all trials
fc = atcm.fun.computefc(s);

% precompute phase-coupling for all trials
[ph,PhsCor] = atcm.fun.computephase(s);

% loop over trial types
for i = 1:length(y)

    % convert firing to Hz/s for this trial
    c(i).fire = f2sr(f{i},DCM.M.dt);
    
    % extract membrane potentials
    c(i).mV   = squeeze(s{i}(:,:,1,:));
    c(i).AMPA = squeeze(s{i}(:,:,2,:)); % AMPA currents
    c(i).GABA = squeeze(s{i}(:,:,3,:)); % GABA-A currents
    c(i).NMDA = squeeze(s{i}(:,:,4,:)); % NMDA currents
    c(i).GABB = squeeze(s{i}(:,:,5,:)); % GABA-B currents
    c(i).Mcur = squeeze(s{i}(:,:,6,:)); % m currents
    c(i).Hcur = squeeze(s{i}(:,:,7,:)); % h currents
    
    c(i).pst = pst;
    c(i).sig = g{i};
    c(i).input = t{i};
    c(i).Spectra = y{i};
    c(i).w       = w;
    
    % layer data is already appropriately stored
    c(i).l = l{i};
    
    % add amplitude (FC) and phase-coupling 
    c(i).phase = ph{i};
    c(i).phase_coupling = PhsCor(i);
    c(i).fc    = fc(i);
    
    
    
end

components = c;