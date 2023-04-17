function [DD,Q] = HVG_DD(x,w,wt)
% returns the Degree Distribiution of the Horizontal Visibility Graph Adjacency
% matrix:
%    [DD,AdjMat] = atcm.fun.HVG_DD(x,w)
% AS22

if nargin < 3 || isempty(wt)
    wt=2;
end

Q  = fast_HVG(x,w,wt);
DD = (1/3) * ( (2/3).^(sum(~~abs(Q))-2) );
DD = DD(:);