function X0 = gauregmat(w)

for i = 1:length(w); 
    X0(i,:) = atcm.fun.makef(w,w(i),4,1);
end

%X0 = X0*X0';