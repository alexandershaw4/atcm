function X = makef(w,Fq,Amp,Wid)
%
% afit.makef(w,Fq,Amp,Wid)
%
% e.g.
% w = 1:100
% S = afit.makef(w,[10 20 50 70],[10 8 5 3],[1 2 5 3]);
% figure;plot(w,S)
%
% AS

if length(Fq) > 1
    for i = 1:length(Fq)
        try
            X0 = X0 + afit.makef(w,Fq(i),Amp(i),Wid(i));
        catch
            X0 =      afit.makef(w,Fq(i),Amp(i),Wid(i));
        end
        %X0(i,:) = afit.makef(w,Fq(i),Amp(i),Wid(i));

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
f   = findthenearest(Fq,w);
f   = f(1);

w   = w - mw;
X   = X + Amp * exp( -(w-f).^2 / (2*(2*Wid)^2) );
w   = w + mw;


