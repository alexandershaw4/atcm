function y = tightenvelope(x,n,w)

dx = real(x);
x  = real(x);
%w  = 1:length(x);
lx = log(x); lx(isinf(lx)) = log(1e-8);
lx=real(lx);
c  = fit(  (w'), (lx),'poly1');
x  = exp( log(x) - c( (w)) );
x = real(x);

x = flipud(filter(ones(1,8)/8, 1, flipud(x))) ;

x = flipud(filter(ones(1,8)/8, 1, flipud(x)));

y=real(x);


%f = fit(w.',x,'gauss4');
%y = f(w);

%[~,I] = findpeaks(x,'NPeaks',n);

% w0 = w(1)-1:w(end)+1;
% for i = 1:length(I)
%     
%     these = [I(i)-1 I(i) I(i)+1];
%     
%     [amp(i),ind0] = max(dx(these));
%     ind(i) = ind0(1);
%      
%     nf(i) = these(ind0);
% end
% 
% y = atcm.fun.makef(w,w(nf),amp,repmat(1.2,length(I),1));
% 
% y=y';