function dx = awinsmooth(x,n)

if nargin < 2 || isempty(n)
    n =1;
end

for i = 1:n
    x  = x(:);
    fx = [x(1); x(1:end-1)];
    bx = [x(2:end); x(end)];
    x  = mean([fx x bx]');
end

dx = x;