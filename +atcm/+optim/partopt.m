function partopt(P,DCM,routine,niter)

% first pass with no reduction
[X0,X,cm,F,j]   = atcm.optim.clusteroptimse(P,DCM,'none',1);
jac             = j{1};
jac(isnan(jac)) = 0;
[~,nc]          = PEig90(real(pinv(jac)));

DCM.M.Nmax = 10;

% Loop
for i = 1:niter
    [X0,X,cm,F,j]   = atcm.optim.clusteroptimse(P,DCM,routine,nc);
    jac             = j{1};
    jac(isnan(jac)) = 0;
    [~,nc]          = PEig90(real(pinv(jac)));
    
    % store 
    X1{i} = X0;
    TM{i} = cm;
    F1(i) = F;
    J1{i} = j;
    
    % update
    P = spm_unvec(X0,P);
    
end