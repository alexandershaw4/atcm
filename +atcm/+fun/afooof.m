function [yfit,win,c,cc] = afooof(w,y,c)

if nargin < 3 || isempty(c)
    % constrained underlying curve fit
    c = atcm.fun.c_oof(w,(y));
end

% Gaussian bump fit
%g = fittype('abs(a1)*exp(-((x-abs(b1))/c1)^2) + abs(a2)*exp(-((x-abs(b2))/c2)^2) + abs(a3)*exp(-((x-abs(b3))/c3)^2)');

g = (y(:)-c(:));

warning off;
win = fit(w.',log(g),'Gauss7');
warning on;

yfit = exp(win(w)) + c;

if nargout > 3

cc(1,:) = win.a1*exp(-((w-win.b1)/win.c1).^2);
cc(2,:) = win.a2*exp(-((w-win.b2)/win.c2).^2);
cc(3,:) = win.a3*exp(-((w-win.b3)/win.c3).^2);



end