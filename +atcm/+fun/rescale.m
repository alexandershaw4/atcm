function y = rescale(x,S1,S2)

if nargin < 3 && length(S1) == 2
    S2 = S1(2);
    S1 = S1(1);
end

y = S1 + (S2-S1) .* (x - min(x) ) / ...
    ( max(x) - min(x) );

end