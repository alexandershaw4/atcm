function y = lsq_widfit(x)


% Gaussian estimation
[y,m] = atcm.fun.lsq_smooth(x);

% Gaussian Parameters
w = (1:length(x))';
f  = m.centres;
opts = m.g(:,m.centres);
b = m.b;

for i = 1:size(opts,2)
    
    ix = opts(:,i);
    on = find(ix);
    
    new(:,i) = interp1(w(on),ix(on),w,[],'extrap');
    new(:,i) = atcm.fun.HighResMeanFilt(new(:,i),1,8);

end

b = atcm.fun.lsqnonneg(new,x(:));
y = new*b;

% % WidFun
% widfun = @(X) sum( (x(:) - sum(smoother(X,opts),2) ).^2 );
% 
% [X,F] = fminsearch(widfun,W);
% 
% end
% 
% function y = smoother(w,opts)
% 
% for i = 1:size(opts,2)
%     x  = opts(:,i);
%     try
%         dx = atcm.fun.HighResMeanFilt(x,1, w(i) );
%     catch
%         dx = x*inf;
%     end
%     y(:,i) = dx;
% end
% 
% end
