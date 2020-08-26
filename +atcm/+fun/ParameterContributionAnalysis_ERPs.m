function yy = ParameterContributionAnalysis_ERPs(DCM,param)
%
% yy = ParameterContributionAnalysis_ERPs(DCM,param)
% 
% pass a fully specified DCM - e.g.
% load FixNew_-101-110_LFPs_PlasticityMean.mat
%
% Select parameter - e.g.
% param = 'CV';


P  = DCM.Ep;
i  = spm_fieldindices(P,param);
pv = spm_vec(P);
f  = @(p) full(spm_cat(atcm.fun.evalERP(DCM,spm_unvec(p,P),DCM.Eg)'));

% perturbation value - relative to posterior (fitted) parameter value
%s  = [.2 .4 .6 .8 1 1.2 1.4 1.6 1.8];
s = [.01 1 100];

% loop and evaluate outputs
for j  = 1:length(i)
    for k = 1:length(s)
        dp           = pv;    
        dp(i(j))     = dp(i(j)) * s(k);     
        yy{j}(k,:)   = f(dp);
    end
end

% name = {'ss' 'sp' 'si' 'dp' 'di' 'tp' 'rt' 'rl'};
% figure;
% for j = 1:8
%     subplot(4,2,j); plot( yy{j}' );
%     title(name{j});
%     if j == 1
%         legend({'.01' '1' '100'});
%     end
% end