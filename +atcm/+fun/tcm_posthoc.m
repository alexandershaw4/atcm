function [p,d] = tcm_posthoc(list,plots)

d.list = list;

for i    = 1:length(list)
    
    load(list{i},'DCM','M');
    
    F(i) = DCM.F;

    fprintf('collecting parameters & evaluating model - %d/%d\n',i,length(list));

    % parameters
    p.Ep(i)     = DCM.Ep;
    p.PP(i,:)   = spm_vec(DCM.Ep);
    p.pE(i)     = DCM.M.pE;
    p.pC(i)     = DCM.M.pC;
    p.V         = diag(spm_vec(DCM.M.pC));
    p.CV(i,:,:) = DCM.Cp;


    % evaluate posteriors and retrieve model prediction and data
    [y,w,s,g,drive,pst,l,oth] = feval(DCM.M.IS,p.Ep(i),DCM.M,DCM.xU);

    % spectrum
    d.dataspec(i,:)  = spm_vec(DCM.xY.y{:});
    d.modelspec(i,:) = spm_vec(y);
    
    % timeseries - squeeze respectively out of state-space
    s0 = squeeze(s{1}(1,:,:,:));

    d.mV(i,:,:)    = squeeze(s0(:,1,:));
    d.ampa(i,:,:)  = squeeze(s0(:,2,:));
    d.gabaa(i,:,:) = squeeze(s0(:,3,:));
    d.nmda(i,:,:)  = squeeze(s0(:,4,:));
    d.gabab(i,:,:) = squeeze(s0(:,5,:));
    d.m_cur(i,:,:) = squeeze(s0(:,6,:));
    d.h_cur(i,:,:) = squeeze(s0(:,7,:));

    % local field potential (e.g. after projecting through J & L)
    d.LFP(i,:) = spm_vec(g);

    % exogenous input
    d.input(i,:) = spm_vec(oth.drive);

    % layer-specific spectra
    d.layerspec(i,:,:) = squeeze(l{1}.weighted);

    % frequencies of psd
    d.w = w(:);

end

if nargin > 2 && (plots)

    figure,
    




end