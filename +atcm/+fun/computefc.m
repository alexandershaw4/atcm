function [FC] = computefc(s)
% takes the integrated model output (states timseries) and computes the
% functional connectivity among states
%
%

for i  = 1:length(s) % loop trials
    s0 = s{i};
    
    % compute correlations between all states in all regions
    [ns,np,nk,nt] = size(s0); 
    dat = reshape(s0, [ns*np*nk nt]);
    
    [R,p]   = corr(dat'); 
    FC(i).R = R;
    FC(i).p = p;
    
end