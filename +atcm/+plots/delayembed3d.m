function Y0 = delayembed3d(data0,delay,col,holdon,nfig)
% performs delay embedding of the inputted chan*time series, using 'delay'
%
%
% input data strcture returned by atcm.fun.iunpack
% ans xv,yv and zv are field names
%
% AS
if nargin < 5 || isempty(nfig)
    nfig = 1;
end

if nargin < 4 || isempty(holdon)
    holdon = 0;
end

thecol = 'w';
if nargin < 3 || isempty(col)
    col = 1;
    thecol = 'w';
elseif ischar(col) 
    thecol = col; 
    col    = 0;
elseif isnumeric(col) && length(col)==3
    thecol = col;
    col = 0;
else
    %col    = 0;
    thecol = 'w';
end

%figure('position',[521         119        1836         866]);
if nfig
    figure('Name','AO','Color',[.3 .3 .3],'InvertHardcopy','off','position',[521 119 1836 866]);
    set(gcf, 'MenuBar', 'none');
    set(gcf, 'ToolBar', 'none');
    drawnow;
end

embed_dimen = 3;

for i = 1:8
    
    delta = delay;
    data  = data0(i,:);
    
    % compute embedding
    data_size=size(data);
    %First we turn data into a row vector.
    if data_size(1)<data_size(2)
        data=data';
    end

    %Then we create a matrix  with circular shift that has our delayed vectors
    %and some extra vectors that we don't need
    data_size=size(data);
    Z=zeros(data_size(1),embed_dimen);


    for j=1:embed_dimen      
        Z(:,j)=circshift(data,[(j-1)*delta,0]);     
    end

    %The extra vectors are erased and the final result is the delays vectors
    %with delay delta and embedding dimension emb_dimen
    Y=Z(delta*embed_dimen-1:end,:);


    % Plot output    
    s(i) = subplot(2,4,i);
    if isnumeric(thecol)
        p = plot3( Y(:,1) , Y(:,2) , Y(:,3)  , 'color', thecol , 'linewidth' , 1 );
    else
        p = plot3( Y(:,1) , Y(:,2) , Y(:,3)  , thecol , 'linewidth' , 1 );
    end
    
    if col == 1
        %f = @jet;
        f = @hot;
        cb = [uint8(f(length(Y))*255) uint8(ones(length(Y),1))]';
        drawnow;set(p.Edge, 'ColorBinding','interpolated', 'ColorData',cb);
    else
    end
    
    Y0{i} = Y;
    
%     xx  = [x(i,:); x(i,:)]; yy = [y(i,:); y(i,:)]; zz = [z(i,:); z(i,:)];
%     cl  = [1:length(x);1:length(x)];
%     hs=surf(xx, yy, zz, cl,'EdgeColor','interp','FaceColor','none','Marker','o') ;
%     shading flat 
    
    %title(c{i});
    %xlabel('time'); ylabel(xv); zlabel(yv); 
    axis square;
    view(3);
    grid on;
    
    s(i).YColor = [1 1 1];
    s(i).XColor = [1 1 1];
    s(i).ZColor = [1 1 1];
    s(i).Color  = [.3 .3 .3];
    
    if holdon
        hold on;
    end
    
end

set(findall(gcf,'-property','FontSize'),'FontSize',20);

