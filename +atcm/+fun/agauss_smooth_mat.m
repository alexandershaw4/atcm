function [P,fm,model] = agauss_smooth_mat(x,n,method,numcomp)
% this function performs sliding-window smoothing with an n-width Gaussian 
% kernel on the input vector x (of length k), stacking each window into a matrix of
% size kxk (a sort of Gaussian-process-like matrix). It then estimates the width of
% each Gauss component and reforms the kxk matrix with appropriate widths.
% Finally it uses non-negative least squares to fit this 'Gaussian matrix'
% back to the original data, eliminating unncessary components. 
% It is a fully linearised Gaussian fitting routine.
%
% [P,matrix] = agauss_smooth_mat(x,n,[m])
%
% 
%
% AS2023

x   = x(:);
w   = (1:length(x))';
fun = @(Wid,f) 2 * exp( -(w-f).^2 / (2*(2*Wid)^2) );

% Gaussian window smoothing
%-----------------------------------------------------------------------
for i = 1:length(x)
    ym(i,:) = ( x.*fun(n,i) );
end

% Estimate widths
%----------------------------------------------------
for i = 1:length(x)
    wd(i) = sqrt( sum(ym(i,:).*(w' - w(i)).^2) /  sum(ym(i,:)) );
end

% Re-estimate matrix with both mu and widths
%-----------------------------------------------------------------------
for i = 1:length(x)
    fm(i,:) = ( ym(i,:)'.*fun(wd(i),i) );
    %fm(i,:) = atcm.fun.makef(w,i,ym(i,i),wd(i));
end

%fm = fm';
   
% Now reduce from a gaussian at every point to only a subset...
%-----------------------------------------------------------------------
if nargin < 3 || method == 1

    % much simpler mixing - non-neg lsq
    b = atcm.fun.lsqnonneg(fm,x);
    i = find(b);
    P = fm.*b';
    P = P(:,i)';   

    % Return components too;
    model.mu     = i(:);
    model.width  = wd(i)';
    
    for k = 1:length(i)
        model.height(k) = P(k,i(k));
    end

elseif method == 2

    % decompose matrix using svd of covariance
    [u,s,v] = svd(cov(fm'));

    i  = atcm.fun.findthenearest(cumsum(diag(s))./sum(diag(s)),.99);
    P  = fm*u(:,1:i);
    P  = abs(P);
    uu = u(:,1:i);

    % estimate mixing for remaining components
    %----------------------------------------------------
    c = x'/abs(P');
    %c = atcm.fun.lsqnonneg(P,x);

    for i = 1:length(c)
        P(:,i) = ( c(i) * P(:,i)' );
    end

    P = P';

elseif method == 3

    % NNMF: works but is NOT deterministic! 
    if nargin < 4 || isempty(numcomp)
        numcomp = 8;
    end

    [W,H] = nnmf(fm,numcomp);

    P = W';

    %c = x'/abs(P);
    c = atcm.fun.lsqnonneg(abs(P)',x);

    for i = 1:length(c)
        P(i,:) = ( c(i) * P(i,:) );
    end

end


% % Now find and eliminate components fully within (under) another component
% %-----------------------------------------------------------------------
% for i = 1:size(P,1)
%     for j = 1:size(P,1)
%          xp = P(i,:);
%          yp = P(j,:);
% 
%          if i~=j && abs(sum(xp - yp)) < 1e-2
%              yp = yp*0;
%              P(j,:) = yp;
%          end
%     end
% end
% 
% r = find(sum(P,2)~=0);
% P = P(r,:);
% 
% % final rescle to componsate for those removed
% %-----------------------------------------------------------------------
% c = atcm.fun.lsqnonneg(abs(P)',x);
% i = find(c);
% P = P(i,:);
% c = c(i);
% 
% for i = 1:length(c)
%     P(i,:) = ( c(i) * P(i,:) );
% end

end









    %wd(i,:) = exp( -(w-i).^2 / (2*(2)^2) )\x;

% b = b(i);
% 
% % optionally reduce to N-components
% %----------------------------------------------------
% if nargin > 2 && ~isempty(m)
%     [~,I] = sort(b,'descend');
%     P = P(I(1:m),:);
%     b = atcm.fun.lsqnonneg(P',x);
%     P = b.*P;
% end


% % decompose matrix using svd of covariance
% %----------------------------------------------------
% [u,s,v] = svd((fm'));
% 
% i  = atcm.fun.findthenearest(cumsum(diag(s))./sum(diag(s)),.99);
% P  = u(:,1:i)'*fm';
% P  = abs(P);
% r  = corr(P',x).^2;
% uu = u(:,1:i);
% 
% % estimate mixing for remaining components
% %----------------------------------------------------
% c = x'/abs(P);
% 
% for i = 1:length(c)
%     P(i,:) = ( c(i) * P(i,:) );
% end
% 
% ym = P;


% repeat?

