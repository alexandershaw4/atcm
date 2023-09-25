function x = replacepeakwithneighbs(x)

if all(size(x) > 1)
    for i = 1:size(x,1)
        x(i,:) = atcm.fun.replacepeakwithneighbs(x(i,:)');
    end
    return;
end

stp = mean(x) + std(x)*3;
stn = mean(x) - std(x)*3;
pk  = [find(x > stp); find(x < stn)];

for i = 1:length(pk)
    
    if pk(i) < length(x) && pk(i) > 1
        x(pk(i)) = 0.5*x(pk(i)-1) + 0.5*x(pk(i)+1);
    elseif pk(i) == 1
        x(pk(i)) = x(pk(i)+1);
    elseif pk(i) == length(x)
        x(pk(i)) = x(pk(i)-1);
    end

end
