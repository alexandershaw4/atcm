function PCSD(X,Y,varargin)
% Plot A/CSD from 3D double input: x(nf,ns,ns)
%            
% Examples:
%           PCSD(X);          % plot 3D double  x [real part]
%           PCSD(X,Y);        % plot 3D doubles x & y 
%           PCSD(X,Y,'imag'); % imaginary part only
%           PCSD(X,Y,'real'); % real part only 
% AS

try Y                                  ; catch Y  = X*0;    end
try s  = eval(['@' lower(varargin{1})]); catch s  = @real;  end
try Hz = varargin{2}                   ; catch Hz = 1:size(X,1); end 
if  isempty(Y); Y = X*0; end

f = @(x,s,m,n)squeeze(s(x(:,m,n)));

p = size(X,2);
b = 0;

for m = 1:p
    for n = 1:p
        if m == n; cl = 'g'; 
        else       cl = 'b';
        end
        b = b + 1;
        subplot(p,p,b);
        plot(Hz,f(X,s,m,n),cl); hold on;
        plot(Hz,f(Y,s,m,n),'r');
        xlim([Hz(1) Hz(end)]);
    end
end
        