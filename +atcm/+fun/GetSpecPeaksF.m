function P = GetSpecPeaksF(w,y)
% Compute spectral peaks
%
% plot( P.w , P.Spec , 'k', P.AlphaF , P.AlphaA , '*' , P.GammaF , P.GammaA , 'o' )
%
% AS

%f = dir('TCM_VS*.mat'); f = {f.name}';

a = [4 13];
b = [13 30];
g = [40 85];

%dw = w(2)-d(1);

for i = 1:size(y,1)
    
    Y = y(i,:);
    
    %load(f{i})
    %Y = DCM.xY.y{1};
    %w = DCM.xY.Hz;
    
    w0 = findthenearest(a(1),w):findthenearest(a(2),w);
    w1 = findthenearest(b(1),w):findthenearest(b(2),w);
    w2 = findthenearest(g(1),w):findthenearest(g(2),w);
    
    
    [ AlphaAmp(i), AlphaFreq(i) ] = GetPks( Y , w , w0 );
    [ BetaAmp(i), BetaFreq(i) ]   = GetPks( Y , w , w1 );
    [ GammaAmp(i), GammaFreq(i) ] = GetPks( Y , w , w2 );
    
    Spec(i,:) = Y;
end

P.AlphaA = AlphaAmp;
P.AlphaF = AlphaFreq;
P.BetaF  = BetaFreq;
P.BetaA  = BetaAmp;
P.GammaA = GammaAmp;
P.GammaF = GammaFreq;
P.Spec   = Spec;
P.w      = w;

end

function [MA,MF] = GetPks(y,w,iw)

[MA,MF] = pksmax(y,w,@real,iw);

if isempty(MA)
    [MA,MF] = pksmax(y,w,@detrend,iw);
end
if isempty(MA)
    MA=nan;
    MF=nan;
end

end

function [MA,MF] = pksmax(y,w,f,iw)

[A,F]   = dpks(y,w,f,iw);
[MA,MF] = truemax(A,F);

end

function [pA,pFq] = dpks(Y,w,f,iw)

[pA,pFq] = findpeaks(f(Y(iw)),w(iw));

end

function [Ma,Mf] = truemax(A,F)

[~,I] = max( A );
Ma    = A(I);
Mf    = F(I);

end