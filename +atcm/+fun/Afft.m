function [Pf,f] = Afft(signal,fs,varargin)
% Return cross spectral density of continuous time signal sampled at fs.
% Include desired frequency vector to return spline interpolated F and Pf.
%
% Use:  [Pf,f] = Afft(signal,fs)           % signal at fs (1/dt)
%       [Pf,f] = Afft(signal,fs,w);        % spline interp'd FoI
%                               w = 1:100; % Freq range of interest
%
% Note - if your frequency range is not integer spaced then your
% effective fs is not 1./dt but dw./dt.   (where dw = 1./(w(2)-w(1) )
%
% AS

if isempty(varargin); varargin{1} = 0; end

if ndims(signal) == 3
    % for channel by sample by trial data
    for i = 1:size(signal,3); 
        [Pf(:,:,:,i),f] = Afft(squeeze(signal(:,:,i)),fs,varargin{1});
    end
    return
end

% calculate cross spectrum using complex fft
Ns = size(signal,1);
Nf = size(signal,2);
Nf = floor(Nf/2) + 1;
Pf = zeros(Nf,Ns,Ns);

L  = size(signal,2);
f  = fs * (0:(L/2))/L;

% just the spectra
if size(signal,1) == 1
    data1  = signal(1,:);
    data1  = fft(data1);
    data1  = (data1/L);
    L2     = floor(L/2);
    Pf     = data1(1:L2+1);
    if varargin{1}
        Hz  = varargin{1};
        warning off
        try
            SPf = interp1(f,data1(1:L2+1),Hz,'spline');
            i   = SPf < 0;
            SPf(i) = -SPf(i);
        catch
            SPf = Hz*0;
        end
        warning on;
    end
else



for s1 = 1:size(signal,1)
    for s2 = 1:size(signal,1)
        
        data1  = signal(s1,:);
        %data1  = pmcov(data1);
        data1  = fft(data1);
        
        data2  = signal(s2,:);
        %data2  = pmcov(data2);
        data2  = fft(data2);
        
        Pxy    = data1 .* conj(data2);
        Pyx    = data2 .* conj(data1);

%         data1  = abs(data1/L);
%         data2  = abs(data2/L);
%         Pxy    = abs(Pxy  /L);
%         Pyx    = abs(Pyx  /L);
        
        data1  = (data1/L);
        data2  = (data2/L);
        Pxy    = (Pxy  /L);
        Pyx    = (Pyx  /L);
        
        L2          = floor(L/2);
        Pf(:,s1,s1) = data1(1:L2+1);
        Pf(:,s2,s2) = data2(1:L2+1);
        Pf(:,s1,s2) = Pxy  (1:L2+1);
        Pf(:,s2,s1) = Pyx  (1:L2+1);
        
        if varargin{1}
            Hz = varargin{1};
            warning off;
            SPf(:,s1,s1) = interp1(f,data1(1:L2+1),Hz,'spline');
            SPf(:,s2,s2) = interp1(f,data2(1:L2+1),Hz,'spline');
            SPf(:,s1,s2) = interp1(f,Pxy  (1:L2+1),Hz,'spline');
            SPf(:,s2,s1) = interp1(f,Pyx  (1:L2+1),Hz,'spline');
            warning on;
            
            i      = SPf < 0;
            SPf(i) = -SPf(i);
        end
    end     
end

end

if varargin{1}
    %fprintf('Returning CSD splined at input frequencies\n');
    Pf = SPf;
    f  = Hz;
end
