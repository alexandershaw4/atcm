function y = agauss_smooth(x,n,a,model)

if nargin < 4 || isempty(model)
    model = 'gauss';
end

x   = real(x(:));
w   = (1:length(x))';
ox  = x;

[x,indx] = atcm.fun.padtimeseries(x);
w   = (1:length(x))';

if nargin < 3 || isempty(a)
    a = 1;
end

if length(a) == 1 && length(n) == 1

    
    switch model
        case 'gauss'; fun = @(Wid,f) a * exp( -(w-f).^2 / (2*(2*Wid)^2) );
        case 'laplace'; fun = @(Wid,f) a *( 1/(sqrt(2)*Wid)*exp(-sqrt(2)*abs(w-f)/Wid) );
    end
    
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