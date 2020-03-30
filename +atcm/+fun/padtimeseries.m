function [y, xi] = padtimeseries(x)
% pad timeseries beginning and end prior filtering
% using reversed timeseries pre and post
% access the original data using y(xi)
% AS 2015

x  = x(:);
xi = size(x);
dx = fliplr(x');
dx = dx(:);

y = [dx; x; dx];

xi = fliplr(xi)+length(x);
xi = xi(1):1:xi(2);

