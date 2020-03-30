function [Data] = PostHocModelSort(f,thing)
% - input a cell array of model filenames and the thing you want to extract
% - defaults to returning spectra and posterior parameters
%
% [] = PostHocModelSort(f)
% [] = PostHocModelSort(f ,'outputs')
% [] = PostHocModelSort(f ,'parameters')
%
%

if nargin < 2 || isempty(thing)
    thing = {'outputs' 'parameters'};
elseif ischar(thing)
    thing = {thing};
end

for i = 1:length(f)
    clear DCM; load(f{i});
    DCM.Ep = EP;
    DCM.Cp = diag(spm_vec(DCM.M.pC));
    
    % OUTPUTS
    if ismember('outputs',thing)
        [y,w,s,g,t,pst,l,n,f0,QD,Sp] = feval(DCM.M.IS,DCM.Ep,DCM.M,DCM.xU);
        
        for j = 1:length(y)
            Data.RSpec (i,j,:,:,:)   = squeeze( DCM.xY.y{j} );
            Data.MSpec (i,j,:,:,:)   = squeeze(        y{j} );
            Data.States(i,j,:,:,:,:) = squeeze(        s{j} );
            Data.Layers(i,j,:,:,:,:) = squeeze(        l{j} );
            Data.firing(i,j,:,:,:,:) = squeeze(        f0{j});
            
            Data.pst                 = pst;
            Data.Hz                  = w;
        end
    end
    
    % PARAMETERS
    if ismember('parameters',thing)
        Data.priors(i)     = DCM.M.pE;
        Data.posteriors(i) = DCM.Ep;
        Data.pricov(i)     = DCM.M.pC;
        Data.poscov(i)     = spm_unvec(diag(DCM.Cp),DCM.M.pC);
        Data.post_vec(i,:) = full( spm_vec( Data.posteriors(i) ));
        Data.pri_vec(i,:)  = full( spm_vec( Data.priors(i) ));
        
        % Also return version multipled by fixed params:
        Data.FullPosteriors(i)  = atcm.fun.unpack_parameters(DCM,DCM.Ep);
        Data.FullPriors(i)      = atcm.fun.unpack_parameters(DCM,DCM.M.pE);
        
        Data.FullPosteriorsVec(i,:) = spm_vec(Data.FullPosteriors(i));
        Data.FullPriorsVec(i,:)     = spm_vec(Data.FullPriors(i));
        
        if i == 1
            % also return an empty parameter structure for reference
            Data.P = spm_unvec( spm_vec(DCM.Ep)*0, DCM.Ep );
        end
    end
end
        
    
        
        
        
        
        
        
        
        