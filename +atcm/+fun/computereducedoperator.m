function V = computereducedoperator(M,fun)
% V = computereducedoperator(M)
% AS2023

if nargin < 2 || isempty(fun); fun = @sum; end

ix = find(fun(M));
V = sparse(1:length(ix),(ix),1,length(ix),size(M,2));
