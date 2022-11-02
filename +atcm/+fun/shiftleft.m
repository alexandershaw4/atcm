function QQ = shiftleft(Q)

QQ = [[Q(2:end,2:end), zeros(size(Q,1)-1,1)];zeros(1,size(Q,2))];