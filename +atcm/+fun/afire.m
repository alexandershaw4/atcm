function F = afire(x,u,v)


%-Format arguments
%--------------------------------------------------------------------------
if nargin < 3, v = 1; end
if nargin < 2, u = 0; end

%-Approximate integral
%--------------------------------------------------------------------------
x    = (x - u)./sqrt(abs(v));
F    = sqrt(1 - exp(-(2/pi)*x.^2))/2;

i    = x < 0;
F(i) = -F(i);
F    = F + 1/2;