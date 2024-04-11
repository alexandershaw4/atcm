function [s,spec,wt] = asinespectrum(w,t,sig,model)
% Spectral estimation by least squares fitting of a sine combination
%
% [~,specrum] = atcm.fun.asinespectrum(w,t,signal)
%            
% where w = frequency vector for output (Freqs of interest)
%       t = sampletimes for signal
%       signal = signal to compute spectrum for
%
% [s,spec] = asinespectrum(w,t,sig,reduce*)
% AS21

if nargin < 4 || isempty(model)
    model = @sin;
end

for i = 1:length(w)
    s(i,:) = model(2*pi*w(i)*(t./1000));
end

if nargin > 2 && ~isempty(sig)
    sig  = sig(:)';
    %spec = sig*s'./(length(sig)./2);
    
    % or use a least-squares fit
    spec = pinv(s*s')*s*sig';
    
    i       = spec < 0;
    spec(i) = -spec(i);
end

% if nargin == 4 && reduce
% 
%     wt = ones(1,size(s,2));
% 
%     for i = 1:10
%         b  = pinv(s*diag(wt)*s')*s*diag(wt)*sig';
%         wt = (b'*s) .* (~~wt);
%         wt = threshold(wt);
%     end
% 
%     spec  = pinv(s*diag(wt)*s')*s*diag(wt)*sig';
% 
%     i       = spec < 0;
%     spec(i) = -spec(i);
% end
% 
% end
% 
% function y = threshold(x)
% 
% [~,I] = sort(x,'descend');
% 
% %w = cumsum(x(I))./sum(x(I));
% 
% %n = atcm.fun.findthenearest(w,0.8);
% 
% ip = find(x(I)>0);
% 
% n = round(length(ip)*.8);
% 
% y             = x*0;
% y(I(ip(1:n))) = x(I(ip(1:n)));
% 
% end