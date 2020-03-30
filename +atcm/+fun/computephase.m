function [ph,PhsCor] = computephase(s)
% takes the integrated model output (states timseries) and computes the
% phase series for the population membrane potentials
%
% also computes the correlation amongst these phases
%

for i = 1:length(s) % loop trials
    for j = 1:size(s{i},1) % loop regions in model
        for k = 1:size(s{i},3) % loop states
            nodedata       = squeeze(s{i}(j,:,k,:));
            %ph{i}(j,:,k,:) = angle( fft(nodedata) );
            ph{i}(j,:,k,:) = angle( hilbert(nodedata) );
        end
    end
end

if nargout == 2
    for i = 1:size(s) % loop trials
        phts = ph{i}; 

        % number of regions, population (per region) and times
        [ns,np,nk,nt] = size(phts); 

        % compute phs correlation between all populations in all regions
        phs   = reshape(phts,[ns*np*nk,nt]);
        [R,p] = corr(phs'); 

        PhsCor(i).R = R;
        PhsCor(i).p = p;

    end
else
    PhsCor = [];
    return;
end