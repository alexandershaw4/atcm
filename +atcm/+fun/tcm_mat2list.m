function out = tcm_mat2list(H,prefix)

C = {'ss' 'sp' 'si' 'dp' 'di' 'tp' 'rt' 'rl'};

if nargin > 1 && ~isempty(prefix)
    n = [prefix '_'];
else
    n = [];
end

if ndims(H) == 2

    list = [];
    for j = 1:8
        for i = 1:8
            if any(H(i,j))
                name = sprintf([n '%s_to_%s'],C{j},C{i});
                out.(name) = H(i,j);
            end
        end
    end
elseif ndims(H) == 3

    list = [];
    for j = 1:8
        for i = 1:8
            if any(squeeze(H(i,j,:)))
                name = sprintf([n '%s_to_%s'],C{j},C{i});
                out.(name) = squeeze(H(i,j,:));
            end
        end
    end





end