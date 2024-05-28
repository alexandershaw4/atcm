function [b,bv,Mu,V,Cov] = agaussreg(X,y)

n = 0;
for i = 1:size(X,2)
    
    if any(X(:,i))

        n = n + 1;

        %dat = X(:,i);
        %mu = mean(dat);

        GX{i} = atcm.fun.VtoGauss(X(:,i));
    
        fun = @(r) sum( (y(:) - atcm.fun.aregress(GX{i},y,'bayes',r(1),r(2))).^2 );
    
        [X0,F] = fminsearch(fun,[1 1]);
    
        [Mu(i,:),Covx] = atcm.fun.aregress(GX{i},y,'bayes',X0(1),X0(2));

        %Mu(i,:) = Mu(i,:) + mu;
    
        V(i,:) = diag(Covx);

        CovStack{n} = Covx;

    else

        Mu(i,:) = zeros(1,length(y));
        V(i,:) = zeros(1,length(y));;
    end

end

b  = pinv(Mu)'*y;
%b  = atcm.fun.aregress(Mu',y,'MAP');

% estimate autovariance
% aV = diag(cov(atcm.fun.VtoGauss(y)));

bv  = pinv(V)'*y;
%bv  = atcm.fun.aregress(V',y,'MAP');

Cov = 0;
for i = 1:length(CovStack)
    Cov = Cov + b(i) * squeeze(CovStack{i});
end
