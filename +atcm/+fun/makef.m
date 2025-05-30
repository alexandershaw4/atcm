function [X,Xk] = makef(w,Fq,Amp,Wid,model)
% Generate a Gaussian/Cauchy/Lapalce/Gamma distribution / mixture, e.g.
%   
%   atcm.fun.makef(w,Fq,Amp,Wid,model)
%
% e.g.
% w = 1:100
% S = afit.makef(w,[10 20 50 70],[10 8 5 3],[1 2 5 3],'gaussian');
% figure;plot(w,S)
%
% AS2019

if nargin < 5 || isempty(model) % 'gaussian or cauchy
    model = 'gaussian';
else
    model = model;
end

if isstruct(Fq)
    Amp = Fq.Amp;
    Wid = Fq.Wid;
    Fq  = Fq.Freq;
end

if length(Fq) > 1
    for i = 1:length(Fq)
        try
            X0 = X0 + atcm.fun.makef(w,Fq(i),Amp(i),Wid(i),model);
        catch
            X0 =      atcm.fun.makef(w,Fq(i),Amp(i),Wid(i),model);
        end
        try
            Xk(i,:) = atcm.fun.makef(w,Fq(i),Amp(i),Wid(i));
        end

    end
    %X0 = max(X0); 
    X  = X0;
    return;
end


try Wid ; catch Wid = 2; end
try Amp ; catch Amp = 2; end

% offset negative values in w
mw  = min(w);
X   = 0*w;
f   = atcm.fun.findthenearest(Fq,w);

try
    f   = f(1);
catch
    f   = Fq(1);
end
f = w(f);

w   = w - mw;
switch model
    case {'Gauss' 'gauss' 'gaussian'}
        X   = X + Amp * exp( -(w-f).^2 / (2*(2*Wid)^2) );
    case 'cauchy'
        X   = X + Amp./(1+((w-f)/Wid).^2);
    case 'laplace'
        X    = X + Amp * ( 1/(sqrt(2)*Wid)*exp(-sqrt(2)*abs(w-f)/Wid) );
    case 'gamma'
        X    = X + Amp * gampdf(w,f);
end
w   = w + mw;

Xk = X;
end


