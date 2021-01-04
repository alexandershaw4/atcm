function y = tightenvelope(x,n,w,P)

dx = real(x);
x  = real(x);
%w  = 1:length(x);
lx = log(x); lx(isinf(lx)) = log(1e-8);
lx=real(lx);
c  = fit(  (w'), (lx),'poly1');
x  = exp( log(x) - c( (w)) );
x = real(x);

%x = flipud(filter(ones(1,8)/8, 1, flipud(x))) ;

%x = flipud(filter(ones(1,8)/8, 1, flipud(x)));

%x = flipud(filter(ones(1,8)/8, 1, flipud(x)));

y=real(x);
% 
% 
% Pk = atcm.fun.GetSpecPeaksF(w,dx');
% 
% X = x'*0;
% 
% X = X + Pk.AlphaA * exp( -(w-Pk.AlphaF).^2 / (2*(2*exp(P.w(1)))^2) );
% X = X + Pk.BetaA * exp( -(w-Pk.BetaF).^2 / (2*(2*exp(P.w(2)))^2) );
% X = X + Pk.GammaA * exp( -(w-Pk.GammaF).^2 / (2*(2*exp(P.w(2)))^2) );
% 
% y = X';
% 

% ab = [findthenearest(w,4):findthenearest(w,14)];
% bb = [findthenearest(w,14):findthenearest(w,35)];
% gb = [findthenearest(w,40):findthenearest(w,80)];
% 
% warning off;
% fa = fit(w(ab).',dx(ab),'gauss2');
% fb = fit(w(bb).',dx(bb),'gauss3');
% fg = fit(w(gb).',dx(gb),'gauss4');
% warning on;
% 
% X = fa(w) + fb(w) + fg(w);
% 
% y = X;

f = fit(w.',x,'gauss5');

% warning off;
% f.c1 = 2*exp(P.w(1));
% f.c2 = 2*exp(P.w(2));
% f.c3 = 2*exp(P.w(3));
% f.c4 = 2*exp(P.w(4));
% warning on;
%

y = f(w);

% [PkV,I] = findpeaks(x,'MinPeakDistance',4);
% 
% [~,Ii]=sort(PkV,'descend');
% 
% if length(I) > n
%     I = I(Ii(1:n));
% end
% 
% % X = x'*0;
% % amp = dx(I);
% % nf = w(I);
% % for i = 1:length(I)
% %     X = X + amp(i) * exp( -(w-nf(i)).^2 / (2*(2*exp(P.w(i)))^2) );
% % end
% % 
% % y = X';
% 
% %amp = dx(I);
% % 
% % %y = atcm.fun.makef(w,w(I),amp,2*exp(P.w(1:length(I))))';
% % 
% % 
% w0 = w(1)-5:w(end)+5;
% for i = 1:length(I)
%     
%     %these = [I(i)-1 I(i) I(i)+1];
%     
%     these = [ I(i)-4:I(i)+4 ];
%     
%     these(these<1)=[];
%     
%     [amp(i),ind0] = max(x(these));
%     ind(i) = ind0(1);
%      
%     nf(i) = these(ind0);
% end
% 
% [~,ord] = sort( w(nf), 'ascend' );
% 
% X = x'*0;
% 
% P.w = ones(1,length(I));
% 
% for i = 1:length(I)
%     X = X + amp(ord(i)) * exp( -(w-w(nf(ord(i)))).^2 / (2*(2*exp(P.w((i))))^2) );
% end
% 
% %y = atcm.fun.makef((1:length(w))./(w(2)-w(1)),w(nf),amp,2*exp(P.w(1:length(I))));
% 
% y=X';