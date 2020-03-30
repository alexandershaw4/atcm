function [bLLE abLLE bLL] = map_Lyapunov(EP,DCM,ix,lab,fields,v1,v2)
%  map_Lyapunov(P,DCM,ix,labs,fields)
%
% Maps the lyapunov exponent of the 3D system of outputs indexed by 'ix' 
% (with labels, 'labs') by varying parameters 'fields' w.r.t each other
% along vectors v1 & v2
%
% e.g. 
% ix  = [2 10 18];
% lab = {'sp mV' 'sp gAMPA' 'sp gGABA-A'}; 
% 
% fields = {'D0(1)' 'D0(2)'};
% v1 = 0.05:0.05:2;
% v2 = 0.05:0.05:2;
%
% atcm.fun.map_Lyapunov(Qp,DCM,[2 10 18],{'sp mV' 'sp gAMPA' 'sp gGABA-A'},{'D0(1)' 'D0(2)'},0.05:0.05:2,0.05:0.05:2)
%
%
% hint on vector multiplier:
%
% k0 = [ fliplr(1./(1:.5:4))*1 ]
% k1 = [ 1*(1:.5:4) ]
% k  = [k0 k1(2:end) ]

% FOR THE 56-STATE THALAMO-CORTICAL MODEL IN atcm.tc_dev:
%                                        L4 L2 L2 L5 L5 L6 TH TH
% states are ordered y(1:8)   = mV     {[ss sp si dp di tp rt rl]}
%                    y(9:16)  = gAMPA  {[ss sp si dp di tp rt rl]}
%                    y(17:24) = gGABAA {[ss sp si dp di tp rt rl]}
%                    y(25:32) = gNMDA  {[ss sp si dp di tp rt rl]}
%                    y(33:40) = gGABAB {[ss sp si dp di tp rt rl]}
%                    y(41:48) = gM     {[ss sp si dp di tp rt rl]}
%                    y(49:56) = gH     {[ss sp si dp di tp rt rl]}

%                                        1  2  3  4  5  6  7  8
%                                        9 10 11 12 13 14 15 16
%                                       17 18 19 20 21 22 23 24
%                                       25 26 27 28 29 30 31 32
%                                       33 34 35 36 37 38 39 40
%                                       41 42 43 44 45 46 47 48
%                                       49 50 51 52 53 54 55 56



% Also perform bifurcation analysis - 
% - select a parameter and move it through a sequence, re-recording le:
% ix  = [2 10 18];
% lab = {'sp mV' 'sp gAMPA' 'sp gGABA-A'}; 
% 
% k0   = 0.05:0.05:2;%[.2 .5 .8  1  1.2 1.5 1.8];
% k1   = 0.05:0.05:2;%[.2 .5 .8  1  1.2 1.5 1.8];

k0=v1;
k1=v2;


for l  = 1:length(k0)
    for j  = 1:length(k1)

        %px        = DCM.Ep;
        px       = EP;
        %px.(x) = log(exp(eval(['px.' fields{1} ]))*k0(j));  % i.e. exp(  log(exp(px.D0(1))*k1) )*60
        %px.(x) = log(exp(eval(['px.' fields{2} ]))*k1(l));
        
        TP1 = log(exp(eval(['px.' fields{1} ]))+k0(j));  % i.e. exp(  log(exp(px.D0(1))*k1) )*60
        TP2 = log(exp(eval(['px.' fields{2} ]))+k1(l));
        
        eval(['px.' fields{1} '= TP1;']);
        eval(['px.' fields{2} '= TP2;']);

        px = spm_vec(px);
        try
            [abLLE(j,l),bLL(j,l,:)] = atcm.fun.phase_lyapunov(EP,DCM.M,DCM.xU,ix,lab,px,0);
            bLLE(j,l) = nanmean( bLL(j,l,:) );
        catch;
             bLLE(j,l) = 0;
             bLL(j,l,:) = 0;
        end
        
    end
end

DCM.Ep = EP;

% Sample 2 extrema points and generate phase space plots
% p(min)
[mapmin,Imin] = min(bLLE(:));
[minx,miny] = ind2sub(size(bLLE),Imin);

px       = DCM.Ep;
%px.D0(1) = eval(['px.' fields{1} ]) + k0(minx);
%px.D0(2) = eval(['px.' fields{2} ]) + k1(miny);
T1 = eval(['px.' fields{1} ]) + k0(minx);
eval(['px.' fields{1} '= T1;']);

T1 = eval(['px.' fields{2} ]) + k1(miny);
eval(['px.' fields{2} '= T1;']);

px       = spm_vec(px);
fx       = @atcm.fun.integrator_pst;
s0       = feval(fx,spm_unvec(px,DCM.Ep),DCM.M,DCM.xU);
s0       = s0(ix,:);

% p(max)
[mapmax,Imax] = max(bLLE(:));
[maxx,maxy] = ind2sub(size(bLLE),Imax);

px       = DCM.Ep;
%px.D0(1) = eval(['px.' fields{1} ]) * k0(maxx);
%px.D0(2) = eval(['px.' fields{2} ]) * k1(maxy);
T1 = eval(['px.' fields{1} ]) + k0(maxx);
eval(['px.' fields{1} '= T1;']);

T1 = eval(['px.' fields{2} ]) + k1(maxy);
eval(['px.' fields{2} '= T1;']);

px       = spm_vec(px);
s1       = feval(fx,spm_unvec(px,DCM.Ep),DCM.M,DCM.xU);
s1       = s1(ix,:);



figure('Name','LYAP','Color',[.3 .3 .3],'InvertHardcopy','off','position',[794 203 1491 997]);
set(gcf, 'MenuBar', 'none');
set(gcf, 'ToolBar', 'none');

subplot(2,3,[1 2 4 5]);
imagesc(bLLE);
axis square;
colormap(hot);colorbar;

ylabel(fields{2},'color','w');
xlabel(fields{1},'color','w');
%ylabel('Cort -> Thal Delay','color','w');
%xlabel('Thal -> Cort Delay','color','w');

DD1 = exp(eval(['DCM.Ep.' fields{1} ]));
DD2 = exp(eval(['DCM.Ep.' fields{2} ]));

CTvec = DD1+k0;
TCvec = DD2+k1;

set(gca,'xtick',1:length(v1),'xticklabels',CTvec);
set(gca,'ytick',1:length(v2),'yticklabels',TCvec);
rotateXLabels(gca,45);

% change color
ticklabels = get(gca,'YTickLabel');
ticklabels_new = cell(size(ticklabels,1),1);
for i = 1:size(ticklabels,1)
    ticklabels_new{i} = ['\color{white} ' ticklabels(i,:)];
end
% set the tick labels
set(gca, 'YTickLabel', ticklabels_new);

ticklabels = get(gca,'XTickLabel');
ticklabels_new = cell(size(ticklabels,1),1);
for i = 1:size(ticklabels,1)
    ticklabels_new{i} = ['\color{white} ' ticklabels(i,:)];
end
% set the tick labels
set(gca, 'XTickLabel', ticklabels_new);

title('Lyapunov Exponent Map','fontsize',20,'color','w')

% symmetric colorbar
Scl = max(abs(bLLE(:)));
caxis([-Scl Scl]);
colormap(cmocean('balance'));

% the phase space: min
%--------------------------------------------------------------
s  = subplot(2,3,3);

Y  = s0';
p  = plot3( Y(:,1) , Y(:,2) , Y(:,3)  , 'r' , 'linewidth' , 1 );
cf = @hot;
cb = [uint8(cf(length(Y))*255) uint8(ones(length(Y),1))]';
drawnow;set(p.Edge, 'ColorBinding','interpolated', 'ColorData',cb);
axis square;
view(3);
grid on;

s(1).YColor = [1 1 1];
s(1).XColor = [1 1 1];
s(1).ZColor = [1 1 1];
s(1).Color  = [.3 .3 .3];
rotate3d();

title('Phase space representation @ min MLE','color','w');

% the phase space: max
%--------------------------------------------------------------
s  = subplot(2,3,6);

Y  = s1';
p  = plot3( Y(:,1) , Y(:,2) , Y(:,3)  , 'r' , 'linewidth' , 1 );
cf = @hot;
cb = [uint8(cf(length(Y))*255) uint8(ones(length(Y),1))]';
drawnow;set(p.Edge, 'ColorBinding','interpolated', 'ColorData',cb);
axis square;
view(3);
grid on;

s(1).YColor = [1 1 1];
s(1).XColor = [1 1 1];
s(1).ZColor = [1 1 1];
s(1).Color  = [.3 .3 .3];
rotate3d();

title('Phase space representation @ max MLE','color','w');

set(findall(gcf,'-property','FontSize'),'FontSize',16)
