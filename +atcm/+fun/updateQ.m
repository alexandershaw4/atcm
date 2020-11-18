function DCM = updateQ(DCM)

y = DCM.xY.y;
w = DCM.xY.Hz;

for i = 1:length(y)
    n      = size(y{i},1);
    m      = size(y{i},2)*size(y{i},3);
    q      = spm_Q(1/2,n,1);
    q      = diag( .25+(w/w(end)) ) .* q;
    q      = q.*hamming(length(w));
    Q{i,i} = kron(speye(m,m),q);
end
DCM.xY.Q  = spm_cat(Q);
DCM.xY.X0 = sparse(size(Q,1),0);
