function [j,ip,list] = example_jaco_for_parameter(DCM,pE,parameter)
% compute & plot the parameter gradients for a specific parameter set
%
%
%

% Get priors (co)variances
pC = DCM.M.pC;
V  = spm_unvec( spm_vec(pC)*0, pC);

% Increments only of specified parameter
V.(parameter) = V.(parameter) + 1;

DCM.M.pC = V;

% Compute parameter Jacobian
[j,ip] = atcm.fun.jaco(DCM,pE);

% Plot the result
figure('position',[1000         424         869         914]);
imagesc(j(ip,:));

% if using intrinsics, get labels
if strcmp(parameter,'H')
    [list,lind] = atcm.fun.intrinsics_names(DCM.M.pE);
else
    % otherwise use an approximate parameter name (numbered)
    list = atcm.fun.DCMVECNAMES(DCM);
    lind = spm_fieldindices(pE,parameter);
    list = list(lind);
end

set(gca,'ytick',1:length(list),'yticklabels',list);
set(findall(gcf,'-property','FontSize'),'FontSize',20)