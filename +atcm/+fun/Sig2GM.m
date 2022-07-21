function [GMfit,GM] = Sig2GM(y,thr)

% "Gaussian filtering"
%--------------------------------------------------------------------------
% fit a signal with a series of gaussians
%
% [GMfit,GM] = atcm.fun.Sig2GM(y)
%
% returns GMfit and the components x (mu), y (density) and w (width) in GM.
%
% 
% Peaks (and dimensionality) are obtained using negative deflections in y 
% and widths are obtained by fitting (fminsearch) back to the original data.
%
% example 1:
%--------------------------------------------------------------------------
% y = [1 1 1 1 1 1 1 1 2 3 2 1 1 1 2 3 4 5 6 7 8 7 6 5 4 3 2 1];
% [GMfit,GM] = atcm.fun.Sig2GM(y);
% plot(1:length(y),y,GM.x,GM.y,'*',1:length(y),GMfit); 
% legend({'input signal' 'peaks' 'GaussFit'},'FontSize',14)
%
% example 2:
%--------------------------------------------------------------------------
% y = [1 1 1 1 1 1 1 1 2 3 2 1 1 1 2 3 4 5 6 7 8 7 6 5 4 3 2 1];
% [GMfit,GM] = atcm.fun.Sig2GM(y);
% % the GM structure also has a function handle, such that:
% w = (1:length(y))';
% GMfit = GM.f(w,GM.x,GM.y,GM.w);
% 
% AS2022

if nargin < 2 || isempty(thr);
    thr = 3;
end

% Identify peaks using diff(x)
y = y(:);
w = (1:length(y))';
d = diff(y);
pk = y*0;
for i = 1:length(y)-2;
    if d(i) > 0 && d(i+1) < 0;
         pk(i+1) = i;
    else pk(i+1) = 0;
    end;
end

pk = find(pk);
wint = (1:length(w))';

D = cdist(pk,pk);
D = mean(D)/round(length(D)/2);
i = find( D < thr );

pk(i)=[];

% remove those marked as too close
pk(pk==0)=[];

% generate a width (component proportion) function
wfun  = @(x) atcm.fun.makef(wint,wint(pk)-1,y(pk),x,'gaussian');
gfun  = @(x) sum( (y(:) - spm_vec(wfun(x(:)))).^2 );
initw = ones(size(pk));

% solve with fminsearch
options.Display = 'off';
[X,F] = fminsearch(gfun,initw,options);

tol = 1;
X(X<tol)=tol;

X(X>8)=4;

Y = wfun(X);
Y = Y(:);

% returns
GMfit = Y;
GM.x = pk(:);
GM.y = y(pk);
GM.w = X(:);
GM.f = @atcm.fun.makef; % e.g. GM.f(w,GM.x,GM.y,GM.w)
GM.pdf = @(w) GM.f(w,GM.x,GM.y,GM.w); % e.g GM.pdf( (1:100)' )