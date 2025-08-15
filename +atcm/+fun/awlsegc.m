function Pf = awlsegc(x,t,w)
% Convert time domain signal to frequency domain using weighted least
% squares method under Gaussian constraint;
%
%     Pf = awlsegc(x,t,w)
%
% AS2024


% Weighted Least Squares Estimator with Gaussian Constraint
for i = 1:length(w)
    F(i,:) = sin(2*pi*w(i)*(t));
end

G  = atcm.fun.VtoRadialGauss(ones(length(w),1),2,[],0); % ~ GP kernel
%W  = pinv(F'*G*F)*F'*G;
%Pf = (x(:)'*W);
Pf = x*(pinv(F)*G);
end