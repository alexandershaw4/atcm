function M = GaussEye(n,wd)

if nargin < 2 || isempty(wd)
    wd = 2;
end

M = zeros(n,n);
w = (1:n)';

for i = 1:n
    
    if length(wd) == n
        x = atcm.fun.makef(w,i,2,wd(i));
    else
        x = atcm.fun.makef(w,i,2,wd);
    end
    M(:,i) = x;
end