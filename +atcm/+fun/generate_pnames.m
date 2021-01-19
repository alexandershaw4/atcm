function n = generate_pnames(P)
% Helper: make vector of parameter names from structure(s)
%

pn = fieldnames(P);
n  = [];

for p = 1:length(pn)
    if  isstruct(P.(pn{p}))
        % recursive
        n = [n generate_pnames(P)];
    else
        % append numbers to fieldnames
        v     = spm_vec(P.(pn{p}));
        for i = 1:length(v)
            t = {sprintf('%s%i',pn{p},i)};
            n = [n t];
        end
    end
end

n = n';

end