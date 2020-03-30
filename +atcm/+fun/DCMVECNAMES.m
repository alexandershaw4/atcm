 function ModelFields = DCMVECNAMES(DCM,varargin)
% for getting the model field names & size from a dcm (or any structure!).
%
% Alex 2014
% update: can handle cell inputs


% Get Model Structure & Field Names
%------------------------------------------
try            P     = DCM.M.pE; 
catch;     P     = DCM.pE; %end
end


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

ModelFields = n';

end