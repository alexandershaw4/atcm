function plotphasevideo(tx)

tx = angle(fft(tx));


for j = 1:size(tx,2)
    x=cos(tx(:,j));
    y=sin(tx(:,j));
    
    % for plots
    s=linspace(0,2*pi,100);
    cx=cos(s);
    cy=sin(s);
    plot(cx,cy);        hold on;
    %plot(x,y,'o',cx,cy)
    scatter(x,y,80,1:length(x),'filled'); hold off;
    axis([-1 1 -1 1])
    axis square
    drawnow
    
end
