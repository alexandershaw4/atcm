function y = GaussKern(n,f,w)

if length(f) > 1
    y = zeros(n,n);
    for i = 1:length(f)
        y = y + atcm.fun.GaussKern(n,f(i),w(i));
    end
    return;
end

x = atcm.fun.makef( (1:n)', f, ones(size(f)), w);
y = x(:)*x(:)';