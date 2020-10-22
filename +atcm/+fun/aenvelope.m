function [e,c,f] = aenvelope(x,n,dolog)
% compute the envelope of x using local maxima, with bias toward the end of
% x. (e.g. for spectra with inherent 1./f power law)
%
% e = atcm.fun.aenvelope(x,n)
%
% AS2020

nx = size(x,1);
w  = (1:length(x))'./length(x);

% Whiten: 1d poly in log-log space
%----------------------------------------------------
warning off;

if nargin < 3 || isempty(dolog)
    dolog = 1;
else
    dolog;
end

if dolog
    % 1D poly model:
    lx = log(x); lx(isinf(lx)) = log(1e-8);
    c  = fit( log(w), (lx),'poly1');
    wx = exp( log(x) - c(log(w)) );
else
    lx = log(x); lx(isinf(lx)) = log(1e-8);
    c  = fit(  (w), (lx),'poly1');
    wx = exp( log(x) - c( (w)) );
end

% Robust linear model:
%c=robustfit(log(w),lx,'welsch');
%wx=exp((c(1)+c(2)*log(w)));
%wx = exp( log(x) - wx );

% Custom  1./f terms model:
% weq = 'a*( 1./x.^-b )';
% startPoints = [1 1];
% c = fit(w,wx,weq,'Start', startPoints);
% wx = wx - c(w);

f  = c;
%c  = exp(c(log(w)));
c = c(w);
warning on;


% Remove bottom 30% frequency content: assuming slow drifts
%----------------------------------------------------
%wx = atcm.fun.bandpassfilter(wx,2*(1./(w(2)-w(1))),[.2 1]*length(x));
%wx = fft( atcm.fun.bandpassfilter(ifft(wx),2*(1./(w(2)-w(1))),[.3 1]*length(x)) );

%wx = (abs(hilbert(wx)));
% wx = wx.*h;

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