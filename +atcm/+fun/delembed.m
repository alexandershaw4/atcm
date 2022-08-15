function y = delembed(x,d)

x      = x(:);
y(:,1) = x;

for i = 2:d
    
    %t = [x(i:end); zeros(i-1,1)];
    
    t = [x(i:end); x(end-1:-1:(end-i+1))];

    y(:,i) = t;

end