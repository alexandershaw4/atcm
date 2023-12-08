function y = agauss_smooth(x,n,a)

x   = real(x(:));
w   = (1:length(x))';
ox  = x;

[x,indx] = atcm.fun.padtimeseries(x);
w   = (1:length(x))';

if nargin < 3 || isempty(a)
    a = 1;
end

if length(a) == 1 && length(n) == 1

    
    fun = @(Wid,f) a * exp( -(w-f).^2 / (2*(2*Wid)^2) );
    
    for i = 1:length(x)
    
        y(i) = mean( x.*fun(n,i) );
    
    end

elseif length(a) == length(x) && length(n) == 1

    fun = @(Wid,f,a) a * exp( -(w-f).^2 / (2*(2*Wid)^2) );
    
    for i = 1:length(x)
    
        y(i) = mean( x.*fun(n,i,a(i)) );
    
    end

elseif length(a) == 1 && length(n) == length(x)
    
    fun = @(Wid,f,a) a * exp( -(w-f).^2 / (2*(2*Wid)^2) );
    
    for i = 1:length(x)
    
        y(i) = mean( x.*fun(n(i),i,a) );
    
    end

end

y = rescale(y(indx),min(x),max(x));