function DCM = updateQ(DCM)

y = DCM.xY.y;
w = linspace(1,2,length(DCM.xY.Hz));

for i = 1:length(y)
    n      = size(y{i},1);
    m      = size(y{i},2)*size(y{i},3);
    q      = spm_Q(1/2,n,1);
    %q      = q.*~eye(length(q));
    
    %q      = diag( w ) + q;
    %q      = q.*hamming(length(w));
    Q{i,i} = kron(speye(m,m),q);
end
DCM.xY.Q  = spm_cat(Q);
DCM.xY.X0 = sparse(size(Q,1),0);
