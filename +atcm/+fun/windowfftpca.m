function [Ppf,su,s] = windowfftpca(pc,fs,w,n,k)
% Computes PSD of signal (sampled at fs) at the frequencies in w, using a
% moving average of n-windows and s.t. the psd being a rank-k
% Gaussian pca over windows.
%
%   Ppf = atcm.fun.windowfftpca(signal,fs,w,n,k)
%
% AS2023


if nargin < 5 || isempty(k)
    k = 8;
end

if nargin < 4 || isempty(n)
    n = 30;
end


% Obtain the power spectrum by window-fft
%[Ppf,~,Pfm] = atcm.fun.AfftSmooth(detrend(pc),fs,w,n);
[Ppf,Pfm] = atcm.fun.tfdecomp((pc),1/fs,w,2,2);

% project eigenvectors over windows onto Gaussian set [V]
Pfm = squeeze(Pfm);
Pfm = real(Pfm);

[u,s,v] = svd((Pfm));

V  = VtoGauss(ones(length(w),1));
ns = size(s,2);
ns = min(ns,size(u,1));
su = u(:,1:ns);

for iw = 1:size(su,2)

 %   B = abs(V*su(:,iw));
    B = su(:,iw);

    %b = atcm.fun.moving_average(B',5)';

    [b,~,Q] = atcm.fun.approxlinfitgaussian(B,[],[],1);
    


    %[uu,ss,vv] = svd(Qc);
    %I = atcm.fun.findthenearest(cumsum(diag(ss))./sum(diag(ss)),.9);
    %q = sum(uu(:,1:I)'*Qc,1);

    % limit number of bumps per vector 
    %nq = min(4,size(Q{1},1));
    %q = Q{1}(1:nq,:);
    %su(:,iw) = sum(q,1);
    
    su(:,iw) = b;
    
end

nkw = 1:k; 
s   = diag(s);
Ppf = abs(su(:,nkw)*s(nkw));

%Ppf = su(:,nkw)*s(nkw,:)*mean(v)';
%Ppf = abs(Ppf);

%s   = diag(s);
%V   = VtoGauss(ones(length(w),1));
%nkw = 1:k; % 1:atcm.fun.findthenearest(cumsum((s))./sum((s)),R);

%Ppf = abs( s(nkw)'*u(:,nkw)' );
%Ppf = atcm.fun.HighResMeanFilt(Ppf,1,2);

end