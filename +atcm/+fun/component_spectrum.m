function [yp,Pf,c] = component_spectrum(w,y,n)

m = min(y);
y = real(y - m);
c  = atcm.fun.c_oof(w,y);
yy = y(:) - c(:);
dy = yy;

for i = 1:n
    [~,I] = max(dy);
    
    Pf(i,:) = atcm.fun.makef(w,w(I)-1,dy(I),1.6,'gaussian');
    
    dy = dy(:) - Pf(i,:)';
    
    coeff(i,:) = [w(I),yy(I)];
end

yp =  m + ( c(:) + spm_vec(sum(Pf,1)) );


%yp = m + c(:) + atcm.fun.findbestdist(w,coeff(:,1),coeff(:,2),repmat(2.6,[n 1]),yy);