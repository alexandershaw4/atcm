function chunk = generate_spaced_windows(fs,t,len)

nw = 1:floor(length(t) ./ (fs*len));
chunk = [1+(len*fs*(nw-1)); (fs*len)+(len*fs*(nw-1))]';