function dx = smoothdct(signal)

nt    = size(signal,2); % number of time points
if nt > 500; bins = floor( linspace(1, nt, 20) );
else       ; bins = floor( linspace(1, nt, 10) );
end

for i = 1:length(bins) / 2
    chunk     = signal(: , bins(i):bins(i+2) );        
    xdct(i,:) = dct(chunk);
end

% Average windows
dx = (mean(xdct,1));