function chunk = gau_comp_fit(Pf0,pct)

if naergin < 2 || isempty(pct)
    pct = .9;
end


GP = atcm.fun.agauss_smooth_mat(Pf0,1);
[u,s,v] = svd(GP);
sn = atcm.fun.findthenearest( cumsum(diag(s))./sum(diag(s)), pct);

for i = 1:sn; chunk(i,:) = sum( u(:,i)*s(i,i)*v(:,i)' ,1); end