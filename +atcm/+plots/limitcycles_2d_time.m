function limitcycles_2d_time(data,xv,yv)
% input data strcture returned by atcm.fun.iunpack
% ans xv,yv and zv are field names
%
%
figure('position',[521         119        1836         866]);

x = data.(xv);
y = data.(yv);
z = data.pst;
c = data.cells;

% check if it needs repping over cells (is 1xn):
s = @(x) size(x,1);
if s(x) ~= 8 ; x = repmat(x,[8,1]); end
if s(y) ~= 8 ; x = repmat(y,[8,1]); end


for i = 1:8
    subplot(2,4,i);
    p = plot3( z , x(i,:) , y(i,:)  , 'r' , 'linewidth' , 2 );
    
    cb = [uint8(jet(length(x))*255) uint8(ones(length(x),1))]';
    drawnow;set(p.Edge, 'ColorBinding','interpolated', 'ColorData',cb);
    
    
%     xx  = [x(i,:); x(i,:)]; yy = [y(i,:); y(i,:)]; zz = [z(i,:); z(i,:)];
%     cl  = [1:length(x);1:length(x)];
%     hs=surf(xx, yy, zz, cl,'EdgeColor','interp','FaceColor','none','Marker','o') ;
%     shading flat 
    
    title(c{i});
    xlabel('time'); ylabel(xv); zlabel(yv); 
    axis square;
    view(3);
    grid on;
end

set(findall(gcf,'-property','FontSize'),'FontSize',20);

