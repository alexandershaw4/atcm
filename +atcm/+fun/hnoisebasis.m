function noise = hnoisebasis(n,p)

w = (1:n)';

y = rescale(hamming(length(w)),1,2).*( w(:).^2 );
y = y./max(y);

noise = p(1)*y + p(2)*flipud(y);

if length(p) > 2
    pp = p(3:end);
    pp = pp(:);
    np = length(pp);
    noise = p(1)*y + p(2)*flipud(y) + atcm.fun.makef(w,pp,ones(np,1),4*ones(np,1));
end

%noise = noise ./ max(noise);