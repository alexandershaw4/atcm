function Pf = awlsegc(x,t,w)
% Convert time domain signal to frequency domain using weighted least
% squares method under Gaussian constraint;
%
%     Pf = awlsegc(x,t,w)
%
% AS2024


% Weighted Least Squares Estimator with Gaussian Constraint
F  = atcm.fun.asinespectrum(w,t);
G  = VtoGauss(ones(length(w),1),4,[],0); % ~ GP kernel
W  = pinv(F'*G*F)*F'*G;
Pf = (x(:)'*W);

end