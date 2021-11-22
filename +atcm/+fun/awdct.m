function c = awdct(x,w)
% Weighted Discrete Cosine Transform (DCT) with weights on the top few
% coefficients specified in w
%
%  c = atcm.fun.awdct(x,w)
%
% AS2021

y = dct(x);
m = dctmtx(length(y));

[XX,ind] = sort(abs(y),'descend');
ind      = ind(1:length(w));

wm      = m*0;
wm(ind) = m(ind);

c     = zeros(length(w),length(x));
for i = 1:length(w)
    c(i,:) =  w(i)*(y(ind(i)) *m(ind(i),:));
end