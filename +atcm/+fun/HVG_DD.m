function [DD,Q] = HVG_DD(x,w)
% returns the Degree Distribiution of the Horizontal Visibility Graph Adjacency
% matrix:
%    [DD,AdjMat] = atcm.fun.HVG_DD(x,w)
% AS22

Q  = fast_HVG(x,w);
DD = (1/3) * ( (2/3).^(sum(Q)-2) );
DD = DD(:);