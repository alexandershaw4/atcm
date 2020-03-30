function [csd,auto,Hz,t0] = timefreq(pst,sig,tw,w)
%
% aspectro.timefreq(pst,signal,tw,w)
%
% pst = sampletimes of values in signal
% signal = signal(s)
% tw = time chunks, e.g. 0:0.25:100
% w = frequencies of interest, e.g. 4:0.25:100
%
% returns:
% csd = time frequency cross spectrum as csd(times,freqs,sig,sig)
% auto = time frequency autospectrum as auto(times,freqs,sig)
% Hz   = corresponding frequency vector
%



% num signals and times
[ns,nt] = size(sig);

% compute sampling frequency
Fs = 1./(pst(2)-pst(1))*1000;

% convert tw to indices
for i = 1:length(tw)
    i0    = findthenearest(tw(i),pst);
    t0(i) = i0(end); 
end

for i = 1:length(t0)-1
    chunk = sig(:,t0(i):t0(i+1));
    [Pf(i,:,:,:),Hz]  = atcm.fun.AfftSmooth(chunk,Fs,w);
end

csd = Pf;

% extract autospectra
for i = 1:ns
    auto(:,:,i) = squeeze(Pf(:,:,i,i));
end

t0 = t0(1:end-1);