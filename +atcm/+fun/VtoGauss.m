function [G,b,GL] = VtoGauss(Q,w,model)
% given a vector that has converted to a smoothed symmetrical matrix,
% perform explicit conversion to Gaussians
%
% e.g. 
%      x = [-10:1:10];           % a vector of data
%      Q = atcm.fun.AGenQn(x,4); % convert to smoothed matrix
%      Q = .5*(Q+Q');            % ensure symmetric
%      G = atcm.fun.QtoGauss(Q); % convert to Gaussian series
%
% AS22

if nargin < 3 || isempty(model)
    model = 'Gauss';
end

if nargin < 2 || isempty(w)
    w = 4;
end

nQ = length(Q);
G  = zeros(nQ,nQ); 
x  = 1:length(Q);

for i = 1:length(Q)
    
    v = Q(i);
    I = i;
    
    if length(w) == length(Q)
        G(:,i) = atcm.fun.makef(x,I-1,v,w(i),model);
    else
        G(:,i) = atcm.fun.makef(x,I-1,v,w,model);
    end

end

if nargout == 2
    % Reduce if requested
    b = atcm.fun.lsqnonneg(G,diag(Q));
else
    b = [];
end

if nargout == 3
    Q  = G;
    A  = Q .* ~eye(length(Q));
    N  = size(A,1);
    GL = speye(N,N) + (A - spdiags(sum(A,2),0,N,N))/4;
end

end