cd('/Users/alexandershaw/Library/CloudStorage/Dropbox/tcm/schizophrenia');

models = dir('aaTCM*.mat'); models = {models.name}';

for i = 1:length(models)

    load(models{i},'DCM','M')

    fprintf('Extracting dataset %d/%d\n',i,length(models));

    data(i,:) = spm_vec(DCM.xY.y);
    model(i,:) = spm_vec(feval(DCM.M.IS,spm_unvec((M.Ep),DCM.M.pE),DCM.M,DCM.xU));

    if i == 1; w = DCM.xY.Hz; end

    Ep(i,:) = real(spm_vec(M.Ep));
    pE(i,:) = spm_vec(DCM.M.pE);

    pC = diag(spm_vec(DCM.M.pC));
    Cp = DCM.Cp;

    F(i) = DCM.F;

   % DCM.M.pE = spm_unvec( real(spm_vec(DCM.M.pE)),DCM.M.pE);
   % DCM.Ep = spm_unvec( real(spm_vec(DCM.Ep)),DCM.Ep);

    DCMs{i} = DCM;

end

g  = contains(models,'patient');
gg = double(g);
gg(gg==0)=-1;

% refit at group level

% Bayesian group inversion using empirical Bayes
%  FORMAT [DCM,PEB,M] = spm_dcm_peb_fit(DCM,M,field)

M = struct;
M.X(:,1) = ones(length(models),1);
M.X(:,2) = gg;

[DCM,PEB,M] = spm_dcm_peb_fit(DCMs',M,'All')

for i = 1:length(models)
    fprintf('Extracting dataset %d/%d\n',i,length(models));
    peb(i,:)=spm_vec(feval(DCM{i}.M.IS,DCMs{i}.Ep,DCM{i}.M,DCM{i}.xU));
    pEp(i,:) = spm_vec(DCM{i}.Ep);
end


% frequentist stats
R = kRandTest(full(pEp(g==1,:)),full(pEp(g==0,:)),0.05,5000);
S = spm_unvec(R.tseries_corr,DCM{1}.Ep)
sig = find(spm_vec(S));
Si = spm_unvec(1:302,S)
nm = generate_pnames(DCM{1}.Ep);

% figure - errorbars for sig niificant params
pn = ceil(length(sig)/5);
for i = 1:length(sig)
    subplot(5,pn,i);
    aplots.dot_with_error(pEp(g==1,sig(i)),pEp(g==0,sig(i)));
    set(gca,'xtick',1:2,'xticklabel',{'SZ' 'CON'});
    title(nm{sig(i)});
end

ri = spm_unvec(V*(1:74)',Si)

% run peb in just sz to model SAP/NS?
%-------------------------------------------
MS = struct;
MS.X(:,1) = ones(length(models(g==1,:)),1);
MS.X(:,2) = score.SANS;
MS.X(:,2) = MS.X(:,2) - mean(MS.X(:,2));

[DCMsz,PEBsz,MSz] = spm_dcm_peb_fit(DCMs(g==1)',MS,'All')

for i = 1:length(DCMsz)
    Epsans(i,:) = spm_vec(DCMsz{i}.Ep);
end

[rv,pv]=corr(real(Epsans*V),score.SANS');
pv = denan(pv);

par_id_cor_sans = find(pv < .05 & pv > 0)

EpsansV = Epsans*V;

active_pars = nm(find(sum(V')))

active_pars(par_id_cor_sans)

figure,
for i = 1:length(par_id_cor_sans)
    subplot(4,4,i);
    s = scatter(EpsansV(:,par_id_cor_sans(i)),score.SANS,90,'filled');
    s.MarkerFaceAlpha=.5;;
    title(active_pars{par_id_cor_sans(i)});
end

% variance explained per datasert
R2 = diag(corr(data',model')).^2;

% get a matrix parameter-space map to reduced space
V = spm_svd(diag(spm_vec(DCM{1}.M.pC)));

% group identifier
g = contains(models,'patient');

% figures of plots
map = alexcmap(56+1);

figure('position',[1102         602        2410         693]);
for i = 1:length(models)
    s1 = subplot(121);
    plot(w,data(i,:),'color',[map(i+1,:) 1/2]); hold on;

    s2 = subplot(122);
    plot(w,model(i,:),'color',[map(i+1,:) 1/2]);hold on;
end

s1.Color = [.4 .4 .4]*2;
s2.Color = [.4 .4 .4]*2;

subplot(121);
plot(w,mean(data,1),'linewidth',4,'color','k');
title('V1 Power Spectrum');
grid on; grid minor;
subplot(122);
plot(w,mean(model,1),'linewidth',4,'color','k');
title('Model Power Spectrum');
grid on; grid minor;
set(findall(gcf,'-property','FontSize'),'FontSize',18)

