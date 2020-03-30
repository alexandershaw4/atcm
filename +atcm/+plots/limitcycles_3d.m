function limitcycles_3d(data,xv,yv,zv)
% input data strcture returned by atcm.fun.iunpack
% and xv,yv and zv are field names
%
% plots the 2-dimensional phase portrait
%
% AS

figure('position',[521         119        1836         866]);

x = data.(xv);
y = data.(yv);
z = data.(zv);
c = data.cells;

% check if it needs repping over cells (is 1xn):
s = @(x) size(x,1);
if s(x) ~= 8 ; x = repmat(x,[8,1]); end
if s(y) ~= 8 ; x = repmat(y,[8,1]); end
if s(z) ~= 8 ; x = repmat(z,[8,1]); end

DoColor = 0;

for i = 1:8
    subplot(2,4,i);
    
    if ~DoColor
        line( x(i,:), y(i,:), z(i,:) , 'color', 'k' );
    else
        p  = plot3( x(i,:), y(i,:), z(i,:)  , 'r' , 'linewidth' , 2 );
        cb = [uint8(jet(length(x))*255) uint8(ones(length(x),1))]';
        drawnow;set(p.Edge, 'ColorBinding','interpolated', 'ColorData',cb);
    end
        
    title(c{i});
    xlabel(xv); ylabel(yv); zlabel(zv);
    axis square;
    view(3);
    grid on;
end

set(findall(gcf,'-property','FontSize'),'FontSize',20);

