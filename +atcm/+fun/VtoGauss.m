function [G,b,GL,wid] = VtoGauss(Q,w,model)
% Approximate conversion of vector x to symmetric feature matrix Q
% containing a Gaussian for each element in input x.
%
% usage: [G,b,GL,wid] = atcm.fun.VtoGauss(x)
%
% AS22

ampwid = 1;

if nargin < 3 || isempty(model)
    model = 'Gauss';
end

if nargin < 2 || isempty(w)
    w = 4;
end

% prerequisites
nQ = length(Q);
G  = zeros(nQ,nQ); 
x  = 1:length(Q);
Q  = real(Q);

% this is a scaled version of Q from which to estimate widths
QQ = rescale(abs(Q),1,(2.58)*2);
QQ = QQ .* sign(Q);

for i = 1:length(Q)
    
    v = Q(i);
    I = i;

    if ampwid
        vv = QQ(i);
        w  = max(real(vv),1);;    
        wid(i) = w;
    end
    
    if length(w) == length(Q)
        G(:,i) = atcm.fun.makef(x,I-1,v,w(i),model);
    else
        G(:,i) = atcm.fun.makef(x,I-1,v,w,model);
    end

end

% symmatric and normalised by matrix norm
G = (G + G')./2;
G = G ./ norm(G);

% now rescaling to min/max of input
G = reshape(rescale(G(:),min(Q),max(Q)),size(G));

if nargout == 2
    % Reduce if requested
    b = atcm.fun.lsqnonneg(G,diag(Q));
else
    b = [];
end

if nargout >= 3
    Q  = G;
    A  = Q .* ~eye(length(Q));
    N  = size(A,1);
    GL = speye(N,N) + (A - spdiags(sum(A,2),0,N,N))/4;
end

end