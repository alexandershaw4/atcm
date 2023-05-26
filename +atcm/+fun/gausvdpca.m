function [y,b,uQ] = gausvdpca(x,N,q)

if ~any(x)
    y  = x;
    b  = 0;
    uQ = 0;
    return;
end

if nargin < 3 || isempty(q)
    q = 20;
end

if nargin < 2 || isempty(N)
    N = 8;
end

if all(size(x)>1)
    for i = 1:size(x,2)
        y(:,i) = atcm.fun.gausvdpca(x(:,i),N,q);
    end
    return;
end


warning off;

Q       = atcm.fun.AGenQn(x,q);
%Q       = (Q + Q')/2;
%Q       = atcm.fun.HighResMeanFilt(Q,1,2);
[u,s,v] = svd(Q);
uQ      = u(:,1:N)'*Q;

b       = lsqminnorm(uQ',x,0,'nowarn');
y       = uQ'*b;

warning on;

end