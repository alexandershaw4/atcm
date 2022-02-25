function [tf,t] = ps2tf(w,y)

if isreal(y)
    phase = zeros(size(y));
else
    phase = imag(y);
end

dt = 1./150;
t  = dt+(0:dt:((.25)-dt));
f  = @(a,f,p) a * sin(2*pi*f*(t-p));

for i = 1:length(w)
    tf(i,:) = f( y(i), w(i), phase(i) );
end

tf = abs(tf);
tf = tf ./ max(spm_vec(tf));