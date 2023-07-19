function y = iterativepeaks(x)

it = true;
n  = 0;
y  = x*0;

while it
    x0   = x;
    n    = n + 1;
    p{n} = atcm.fun.indicesofpeaks(x);
    y(p{n}) = n;
    x(p{n}) = 0;

    if all(x==x0)
        it = false;
    end
end
