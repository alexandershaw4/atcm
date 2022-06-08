function [y,hx] = tfdecomp(x,dt,w,space,order,method)
% Spectral estimation by time-frequency decomposition of a signal x.
%
% usage: [y,hx] = tfdecomp(x,dt,w,space,order)
%
% AS22

if nargin < 6 || isempty(method); method = @mean; else; method = method; end
if nargin < 5 || isempty(order);  FilterOrder = 4; else; FilterOrder = order; end
if nargin < 4 || isempty(space);  space = 8; else; space = space; end

% generate frequency windows
win   = [w(:)-space/2 w(:)+space/2];
y     = w(:)*0;

% move edges outside range(w) back inside
win(min(win,w(1))~=w(1))    = w(1);
win(max(win,w(end))>w(end)) = w(end);

% symmetrically pad the timeseries such that xpad(ind) = x
[xpad,ind] = atcm.fun.padtimeseries(x);

for i = 1:length(w)
    
    % rgenerate filter coefficients
    [B, A] = atcm.fun.mk_filter(1./dt,win(i,1),win(i,2),FilterOrder);
    
    % pass through filter forward and backward and average
    fx = filtfilt(B,A,xpad(:));
    
    % get hilbert envelope
    hx(i,:) = hilbert(abs(fx));
    
    % extract amplitude for this frequency step
    y(i) = method(hx(i,ind));
    
end

% remove padding
hx = hx(:,ind);


% C = cov(hx);
% [RHO,LAMBDA] = eig(C);
% LAMBDA = diag(LAMBDA);               % extract the diagonal elements
% [LAMBDA,ind]=sort(LAMBDA,'descend'); % sort eigenvalues
% RHO = RHO(:,ind);                    % and eigenvectors
% 
% y = hx*RHO(:,1);


% [u,s,v] = svd(hx');
% 
% yy = u(:,1)'*hx';
% i = yy<0;
% yy(i) = -yy(i);
% y=yy(:);



% for i = 1:size(hx,2)
%     this = hx(:,i);
%     this = atcm.fun.HighResMeanFilt(this,1,4);
%     spec(:,i) = this;
% end
% y = mean(spec,2);

% % apply smoothing to the TF matrix and extract metric
% f  = atcm.fun.HighResMeanFilt(hx,1,4);
% b  = atcm.fun.HighResMeanFilt(flipud(fliplr(hx)),1,4);
% l  = (f + flipud(fliplr(b)))./2;
% hx = hx(:,ind);
% 
% for i = 1:length(w)
%     
%     y(i) = method(l(i,:))./3;
%     
% end






% % or Eigendecomposition
% [u,s,v] = svd(atcm.fun.HighResMeanFilt(hx,1,4));
% s = diag(s);
% i = atcm.fun.findthenearest(cumsum(s)./sum(s),.99);
% s = diag(s);
% l = u(:,1:i)*s(1:i,1:i)*v(:,1:i)';
% 
% for i = 1:length(w)
%     
%     y(i) = method(l(i,:));
%     
% end