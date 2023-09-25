function [comps,I,Q,ve] = approxlinfitgaussian(x,thresh,method,rankg)
% Approximate Linearised Gaussian Fitting to x
%
% usage: [comps,I,Q,ve] = atcm.fun.approxlinfitgaussian(x,thresh)
%
% notes: ordinarily fitting a (set of) Gaussian(s) to some data is a
% nonlinear optimisation that requires fitting the mean (x-position), height 
% (y-position) and width of each Gauss. However under very specific
% assumptions this can be approximated by a linearisation. 
% % Here a few options are available:
%
% Option 1: (default) or set method = 'none'
%
% 1. conver vector y to symmetric matrix Q by placing gaussian in each column
% [i] with mean of x[i], height y[i] and scaled width
% 2. find indices (pk) of 'peaks' in y based on sign -1 gradients
% 3. take dy = max Q(pk,:)
% 4. explicitly convert to a series of Gaussians using obtained means,
% heights and widths.
%
% Option 2: set method = 'corr'
%
%  1. for every poisition in vector y: x[i] produce a Gaussian with mean
%  x[i], amplitude y(i) and width scaled to 2 s.d. of the amplitude
%  2. from thematrix generated in (1) correlate every vector with y to get
%  variance explained (r.^2).
%  3. rank-sort and compute cumulative r.^2 and find thresholh point of
%  n-comps required to explain k-% variance (2nd input to function)
%  4. retain only these vectors of the matrix, so that matrix Q goes from
%  size Q(length(y),length(y)) to size Q(num_comps,length(y))
%  5. from this new component matrix (where each feature vector contains a
%  Gaussian), find which vectors are 'visible' i.e. max, so that Q further
%  reduces by removing features hidden under other features.
%
% Options 3: set method = 'lsq'
%
% 1. least squares fit to Gaussian matrix (e.g. Gauss process)
%
%
% AS2023

if nargin < 2 || isempty(thresh)
    thresh = 0.9;
end
if nargin < 3 || isempty(method)
    method = 'none';
end

% iterative: refit to n-th residuals
if nargin > 3 && ~isempty(rankg)
    xx = x;
    for i = 1:rankg
        [dx(i,:),I{i},Q{i}] = atcm.fun.approxlinfitgaussian(xx,thresh,method);
        xx = xx(:) - dx(i,:)';
    end
    comps = sum(dx,1);
    return;
end


% algorithm
x  = real(x);
[Q,~,~,W]  = atcm.fun.VtoGauss(x);

m = diag(Q) - x;
for i = 1:length(Q)
    Q(i,:) = Q(i,:) - m(i);
end


switch method
    case 'none';
        
        x   = real(x(:));
        [Q,~,~,W]  = atcm.fun.VtoGauss(x);
        pk  = atcm.fun.indicesofpeaks(real(x));
        if length(pk) > 1
            x   = max(Q(pk,:));
        else
            x = Q(pk,:);
        end
        
        w   = 1:length(x);

        if length(W) ~= length(pk)
            W = repmat(W(1),1,length(w));
        end

        % now we have mean, height and width for each gaussian - 
        [comps,Q] = atcm.fun.makef(w,w(pk)-1,x(pk),W(pk)/2);

        I.mu  = w(pk)-1;
        I.amp = x(pk);
        I.wid = W(pk)./2;

    case 'corr'
        ve = corr(x,Q).^2;
        
        [~,I]=sort(ve,'descend');
        
        n = atcm.fun.findthenearest(cumsum(ve(I))./sum(ve),thresh);
        these = I(1:n);
        
        comps = Q(these,:);
        
        % work backwards form this model
        [~,bwi]=max(comps);
        
        comps = comps(unique(bwi),:);



    case 'lsq'

        b = Q\x;
        bb = b;
        bb(bb<0)=0;

        bQ = bb'*Q;
        bQ = bQ ./ max(bb);

        C = bb'.*Q;

        removed = find(var(C)==0);

        C(:,removed) = [];

        comps = C';

        comps = max(x)*(comps./max(comps(:)));

end
    
end