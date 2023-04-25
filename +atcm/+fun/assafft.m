function Y = assafft(x,dt,w,N)
% returns the N-component power spectrum of x at frequencies in w:
%
% Y = atcm.fun.assafft(x,dt,w,N)
%
% using tf-filtering of a Singular Spectrum Analysis of x
% 
% AS2023

%X = atcm.fun.assa(x,N);
X = x;

Pf = atcm.fun.tfdecomp(X,dt,w,8,2);

    %Y(:,i) = atcm.fun.gaufun.SearchGaussPCA(abs(Pf),8);
    
    Q = atcm.fun.AGenQn(Pf,8);
    %Q = Q .* ~eye(length(Q));
    [u,s,v] = svd(Q);
    Y = abs( Q*v(:,1:N) );

% for i = 1;%:N
%     Pf = atcm.fun.tfdecomp(X(:,i),dt,w,8,1);
% 
%     %Y(:,i) = atcm.fun.gaufun.SearchGaussPCA(abs(Pf),8);
%     
%     Q = atcm.fun.AGenQn(Pf,8);
%     [u,s,v] = svd(Q);
%     Y(:,i) = abs( Q*v(:,1:N) );
% 
% end