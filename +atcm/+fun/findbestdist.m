function [opty,m] = findbestdist(w,f,a,wid,y,elimination)
% Given observation x-vector w and positions of peaks on x-axis (f) and 
% corresponding y-axis values (a) and widths (wid), generate and select 
% from a combination of Gaussians, Cauchy, Laplace and Gamma distributions 
% that best represent data in vector y.
%
% opty = atcm.fun.findbestdist(w,f,a,wid,y)
%
% e.g. for a power spectrum y, at frequencies w, specify a peak at 50 hz
% with amplitude 5 and wid 4:
%   atcm.fun.findbestdist(4:90,50,5,4,y)
% the function will fit the distbribution type that explains the most
% variance in input power spectrum y. You can specify any number of peaks 
% and the sum of the best distribution for each will be returned.
%
% AS22

if nargin < 6 || isempty(elimination)
    elimination = 0;
end

if isstruct(f)
    a   = f.Amp;
    wid = f.Wid;
    f   = f.Freq;
end

% generate distbriutions
for i = 1:length(f)
    m0(i,:) = atcm.fun.makef(w,f(i),a(i),wid(i),'gaussian');
    m1(i,:) = atcm.fun.makef(w,f(i),a(i),wid(i),'cauchy');
    m2(i,:) = atcm.fun.makef(w,f(i),a(i),wid(i),'laplace');
    m3(i,:) = atcm.fun.makef(w,f(i),a(i),wid(i),'gamma');
end

mm(1,:,:) = m0;
mm(2,:,:) = m1;
mm(3,:,:) = m2;
mm(4,:,:) = m3;

f  = @(x,yy) lmod(x,mm,yy);
x0 = ones(size(mm,2),1); 

% (hierarchically) find best fits for each component
nc  = size(mm,2);
xx0 = x0;
for i = 1:nc
    for j = 1:4
        if i == 1 || ~elimination
            yy = y;
        else
            yy = yy - squeeze(mm(xx0(i-1),i-1,:));
        end
        
        xx0(i) = j;
        e(j) = f(xx0,yy);
        
        [~,I] = min(e);
        xx0(i) = I;
    end
end
        
for i = 1:nc
    m(:,i) = squeeze(mm(xx0(i),i,:));
end

opty = sum(m,2);

end

function g = lmod(x,mm,y)

    x = (round(x));
   
    mm = real(mm);
    y  = real(y);
    
    for i = 1:length(x)
        m(:,i) = squeeze( mm(x(i),i,:) );

        %e = y(:) - sum(m,2);
        %g = sum( e.^2 );

        g = 1 - (corr(y(:),sum(m,2))).^2;
        
    end
end