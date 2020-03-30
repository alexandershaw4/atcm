function [list,ind] = intrinsics_names(pE)
% returns a vector list of intrinsic connections in the order of
% spm_vec(P.H)
%
%
%
%

C = {'ss' 'sp' 'si' 'dp' 'di' 'tp' 'rt' 'rl'};

n = 0;
for i = 1:length(C)
    for j = 1:length(C)
        n = n + 1;
        list{n} = sprintf('%s -> %s',C{i},C{j});
    end
end
list=list';

% if a full param strcuture is supplied, find indices of H in full list
if nargin > 0 && isstruct(pE) 
    ind = spm_fieldindices(pE,'H');
end