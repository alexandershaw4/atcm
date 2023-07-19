function y = asmooth(x,n)

if nargin > 1
    for i = 1:n
        x = atcm.fun.asmooth(x);
    end
    y = x;
    return;
end

x  = real(x);
x  = x(:);
x0 = [x(2:end);x(end)]./2;
x1 = [x(1);x(1:end-1)]./2;

y = (x+x0+x1)./3;