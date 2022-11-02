function L = symfluctu(w,a,x)

if nargin < 3 || isempty(x)
    x = 4;
end

n = length(a);
b = linspace(w(1),w(end),n+1);
b = b(1:n) + diff(b)/2;
L = atcm.fun.makef(w,b-w(1),a,~~b*x);
