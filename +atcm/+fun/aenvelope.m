function e = aenvelope(x,n)
% compute the envelope of x using local maxima, with bias toward the end of
% x. (e.g. for spectra with inherent 1./f power law)
%
%
% AS2020

nx = size(x,1);
w  = (1:length(x))'./length(x);

% Whiten: 1d poly in log-log space
%----------------------------------------------------
warning off;
lx = log(x); lx(isinf(lx))=log(1e-8);
c  = fit(log(w), (lx),'poly1');
wx = exp( log(x) - c(log(w)) );
warning on;

% Find local maxima
%----------------------------------------------------
if nx > n+1
    %[~,iPk] = atcm.fun.maxpoints(double(x.*w),n,'max');
    [~,iPk] = atcm.fun.maxpoints(double(wx),n,'max');
    iPk = sort(iPk,'ascend')'; 
else
    iPk = [];
end

% Ensure first and last points included
%----------------------------------------------------
iLocs = [1; iPk; nx];
iLocs = unique(iLocs);

% smoothly connect the maxima via a spline
%----------------------------------------------------
e = interp1(iLocs,x(iLocs),(1:nx)','spline'); % 'spline'
e = abs(e);