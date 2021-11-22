function [tf,t] = ps2tf(w,y)

dt = 1./150;
t  = dt+(0:dt:((.25)-dt));
f  = @(a,f) a * sin(2*pi*f*t);

for i = 1:length(w)
    tf(i,:) = f( y(i), w(i) );
end

tf = abs(tf);
tf = tf ./ max(spm_vec(tf));