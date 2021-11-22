function [yt,ct] = spec2series(w,y,t)

t  = t ./ 1000;
dt = t(2) - t(1);
yt = t*0;
t  = t - t(1);

for i  = 1:length(w)
    xt = y(i) * sin(2*pi*w(i)*t + imag(i) );
    yt = yt + [xt(i:end); zeros((i-1),1)];
    ct(i,:) = xt;
end