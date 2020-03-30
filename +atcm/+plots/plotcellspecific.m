function plotcellspecific(units,vect)
% convenient plotting function for the 8 cells - e.g. for time-series of
% spectra or firing ... give unit vector in first input and 8-by-n in
% second
%

cells = {'ss'  'sp'  'si'  'dp'  'di'  'tp'  'rt'  'rl'};

for i = 1:8
    subplot(2,4,i); plot(units,vect(i,:),'linewidth',2);
    title(cells{i});
end