function limitcycles_2d(data,xv,yv,col)
% input data strcture returned by atcm.fun.iunpack
% and xv,yv and are field names
%
% plots the 2-dimensional phase portrait
%
% AS

if nargin < 4 || isempty(col)
    col = 'k';
end

%figure('position',[521         119        1836         866]);

x = data.(xv);
y = data.(yv);
c = data.cells;

% check if it needs repping over cells (is 1xn):
s = @(x) size(x,1);
if s(x) ~= 8 ; x = repmat(x,[8,1]); end
if s(y) ~= 8 ; x = repmat(y,[8,1]); end

for i = 1:8
    subplot(2,4,i);
    line( x(i,:), y(i,:), 'color', col );hold on;
    
%     xx  = [x(i,:); x(i,:)]; yy = [y(i,:); y(i,:)]; zz = [z(i,:); z(i,:)];
%     cl  = [1:length(x);1:length(x)];
%     hs=surf(xx, yy, zz, cl,'EdgeColor','interp','FaceColor','none','Marker','o') ;
%     shading flat 
    axis square;
    title(c{i});
    xlabel(xv); ylabel(yv);
    grid on;
end

set(findall(gcf,'-property','FontSize'),'FontSize',20);

