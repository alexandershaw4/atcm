function p = indicesofpeaks(x)

x = x(:);
d = gradient(x);
s = sign(d);
p = find(s==-1);

% determine beginning of peak / negative descent
for i = 2:length(s)
    
    if (ismember(i,p) && s(i-1) ~= -1) && (ismember(i,p) && s(i-1) ~=0)
        s(i) = -1;
    elseif (ismember(i,p) && s(i-1) == -1) || (ismember(i,p) && s(i-1) == 0)
        s(i) = 0;
    else
        s(i) = 1;
    end

end

p = find(s==-1);

% check we haven't lnded 1-point beyond max
for i = 1:length(p)
    if p(i) > 1
        if x(p(i)) < x(p(i)-1)
            p(i) = p(i) - 1;
        end
    end
end


end