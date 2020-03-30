function [Pf,f,Pfm] = AfftSmooth(signal,fs,varargin)
% A smoothing wrapper on Afft.m for CSDs.
% 
% 
%
%

nt    = size(signal,2); % number of time points
if nt > 500; bins = floor( linspace(1, nt, 20) );
else       ; bins = floor( linspace(1, nt, 10) );
end

try varargin{2} ; 
    bins = floor( linspace(1, nt, varargin{2}) );
end

for i = 1:length(bins) / 2
    chunk  = signal(: , bins(i):bins(i+2) );    
    if nargin > 2
        [Pfm(:,:,:,i),f] = atcm.fun.Affti(chunk,fs,varargin{1});
    else
        w = linspace(1 , size(signal,2)/2 , size(signal,2) );
        [Pfm(:,:,:,i),f] = atcm.fun.Affti(chunk,fs,w);
    end
end

% Average windows
Pf = (mean(Pfm,4));